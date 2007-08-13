/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <sstream>
#include <cassert>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "BGAutils.hpp"
#include "BGAlib.hpp"
#include "AtomSequence.hpp"
#include "AtomCost.hpp"
#include "Counter.hpp"
#include "Random.hpp"
#include "R3linalg.hpp"

RegisterSVNId BGAlib_cpp_id("$Id$");

using namespace LIGA;

// constants

double BGA::eps_badness = 1.0e-10;

////////////////////////////////////////////////////////////////////////
// common helper functions
////////////////////////////////////////////////////////////////////////

// read lines that do not start with number
bool read_header(istream& fid, string& header)
{
    double x;
    string line;
    istringstream istrs;
    header.clear();
    for (   istream::pos_type p = fid.tellg();
	    !fid.eof() && getline(fid, line);
	    p = fid.tellg()
	)
    {
	istrs.clear();
	istrs.str(line);
	if (istrs >> x)
	{
	    fid.seekg(p);
	    break;
	}
	else
	{
	    header.append(line + '\n');
	}
    }
    return !(fid.rdstate() & ios::badbit);
}

bool read_header(istream& fid)
{
    string dummy;
    return read_header(fid, dummy);
}

// read as many numbers as possible
template<typename T> bool read_data(istream& fid, vector<T>& v)
{
    // prepare v
    T x;
    while (fid >> x)
    {
	v.push_back(x);
    }
    return !(fid.rdstate() & ios::badbit);
}


double penalty(double dd)
{
    static Counter* penalty_calls = Counter::getCounter("penalty_calls");
    penalty_calls->count();
    return dd*dd;
}

////////////////////////////////////////////////////////////////////////
// Atom_t definitions
////////////////////////////////////////////////////////////////////////

Atom_t::Atom_t(double rx0, double ry0, double rz0, double bad0) :
    fixed(false), ttp(LINEAR), badness(bad0)
{
    r[0] = rx0;
    r[1] = ry0;
    r[2] = rz0;
}

double Atom_t::Badness() const
{
    return badness;
}

double Atom_t::FreeBadness() const
{
    return fixed ? 0.0 : badness;
}

double Atom_t::IncBadness(double db)
{
    badness += db;
    return badness;
}

double Atom_t::DecBadness(double db)
{
    badness -= db;
    if (badness < 0.0)	badness = 0.0;
    return badness;
}

double Atom_t::ResetBadness(double b)
{
    badness = b;
    return badness;
}

bool operator==(const Atom_t& a1, const Atom_t& a2)
{
    return equal(a1.r.data(), a1.r.data() + 3, a2.r.data());
}

double dist2(const Atom_t& a1, const Atom_t& a2)
{
    static Counter* distance_calls = Counter::getCounter("distance_calls");
    distance_calls->count();
    double dr[3] = { a1.r[0]-a2.r[0], a1.r[1]-a2.r[1], a1.r[2]-a2.r[2] };
    return dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
}


////////////////////////////////////////////////////////////////////////
// AtomFilter_t - definitions of subclasses
////////////////////////////////////////////////////////////////////////

bool BondAngleFilter_t::Check(Atom_t* pta, Molecule* pm)
{
    // first neighbors are closer than max_blen
    list<Atom_t*> first_neighbors;
    list<double>  first_distances;
    // second neighbors are closer than 2*max_blen
    list<Atom_t*> second_neighbors;
    list<double>  second_distances;
    // find first and second neighbors and corresponding distances
    typedef vector<Atom_t*>::iterator VPAit;
    for (AtomSequence seq(pm); !seq.finished(); seq.next())
    {
	double blen = dist(*pta, *seq.ptr());
	if (0.0 < blen && blen < 2*max_blen)
	{
	    second_neighbors.push_back(seq.ptr());
	    second_distances.push_back(blen);
	    if (blen < max_blen)
	    {
		first_neighbors.push_back(seq.ptr());
		first_distances.push_back(blen);
	    }
	}
    }
    // check all new bond angles formed by test atom
    typedef list<Atom_t*>::iterator LPAit;
    typedef list<double>::iterator LDit;
    LPAit pfn1i, pfn2i;
    LDit dfn1i, dfn2i;
    for (   pfn1i = first_neighbors.begin(), dfn1i = first_distances.begin();
	    pfn1i != first_neighbors.end();  ++pfn1i, ++dfn1i )
    {
	// check test atom angles, ta_angle
	pfn2i = pfn1i; ++pfn2i;
	dfn2i = dfn1i; ++dfn2i;
	for (; pfn2i != first_neighbors.end(); ++pfn2i, ++dfn2i)
	{
	    double dsqneib = dist2(**pfn1i, **pfn2i);
	    double ta_angle = 180.0 / M_PI * acos(
		    ( (*dfn1i)*(*dfn1i) + (*dfn2i)*(*dfn2i) - dsqneib ) /
		    (2.0*(*dfn1i)*(*dfn2i))  );
	    if (ta_angle < lo_bangle || ta_angle > hi_bangle)
	    {
		return false;
	    }
	}
	// check neighbor atom angles, na_angle
	LPAit psni = second_neighbors.begin();
	LDit dsni = second_distances.begin();
	for (; psni != second_neighbors.end(); ++psni, ++dsni)
	{
	    double dneib = dist(**pfn1i, **psni);
	    if (dneib == 0.0 || !(dneib < max_blen))
	    {
		continue;
	    }
	    double na_angle = 180.0 / M_PI * acos(
		    ( (*dfn1i)*(*dfn1i) + dneib*dneib - (*dsni)*(*dsni) ) /
		    (2.0*(*dfn1i)*dneib)  );
	    if (na_angle < lo_bangle || na_angle > hi_bangle)
	    {
		return false;
	    }
	}
    }
    return true;
}

bool LoneAtomFilter_t::Check(Atom_t* pta, Molecule* pm)
{
    // atom is always good with respect to empty molecule
    if (pm->NAtoms() == 0)
    {
	return true;
    }
    // find whether any atom in the molecule is closer than max_dist
    bool has_buddy = false;
    for (AtomSequence seq(pm); !seq.finished() && !has_buddy; seq.next())
    {
	double d = dist(*pta, *seq.ptr());
	has_buddy = (0.0 < d) && (d < max_dist);
    }
    return has_buddy;
}


////////////////////////////////////////////////////////////////////////
// DistanceTable definitions
////////////////////////////////////////////////////////////////////////

DistanceTable::DistanceTable() : vector<double>()
{
    init();
}

DistanceTable::DistanceTable(const double* v, size_t s) : vector<double>()
{
    resize(s);
    copy(v, v+s, begin());
    init();
}

DistanceTable::DistanceTable(const char* file) : vector<double>()
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    bool result = read_header(fid) && read_data(fid, *this);
    // check if everything was read
    if ( !result || !fid.eof() )
    {
	fid.clear();
	cerr << "E: " << file << ':' << fid.tellg() <<
	    ": error reading DistanceTable" << endl;
	throw IOError();
    }
    fid.close();
    init();
}

DistanceTable::DistanceTable(const vector<double>& v) : vector<double>(v)
{
    init();
}

DistanceTable::DistanceTable(const DistanceTable& d0)
{
    *this = d0;
}

DistanceTable& DistanceTable::operator= (const vector<double>& v)
{
    vector<double>& this_vector = *this;
    this_vector = v;
    init();
    return *this;
}

DistanceTable& DistanceTable::operator= (const DistanceTable& d0)
{
    if (this == &d0)
	return *this;
    resize(d0.size());
    copy(d0.begin(), d0.end(), begin());
    NAtoms = d0.NAtoms;
    Nuniqd = d0.Nuniqd;
    max_d = d0.max_d;
    return *this;
}

DistanceTable::iterator DistanceTable::find_nearest(const double& dfind)
{
    iterator ii = lower_bound(begin(), end(), dfind);
    if (    ( ii == end() && size() != 0 ) ||
	    ( ii != begin() && (dfind - *(ii-1)) < (*ii - dfind) )
       )
	--ii;
    return ii;
}

DistanceTable::iterator
DistanceTable::find_nearest_unused(const double& d, valarray<bool>& used)
{
    iterator ii = find_nearest(d);
    int idx = ii - begin();
    if (used[idx])
    {
        int hi, lo, nidx = -1;
        int sz = size();
        for (hi = idx + 1; hi < sz && used[hi]; ++hi)
        { }
        if (hi < sz)
        {
            nidx = hi;
        }
        for (lo = idx - 1; lo >= 0 && used[lo]; --lo)
        { }
        if (lo >= 0 && (nidx < 0 || d - at(lo) < at(nidx) - d))
        {
            nidx = lo;
        }
        idx = nidx;
        ii = begin() + nidx;
    }
    return ii;
}

DistanceTable::iterator DistanceTable::return_back(const double& dback)
{
    iterator ii = lower_bound(begin(), end(), dback);
    return insert(ii, dback);
}

vector<double> DistanceTable::unique()
{
    vector<double> dtu(Nuniqd);
    vector<double>::iterator dui = dtu.begin();
    double eps_dd = sqrt(BGA::eps_badness);
    double d0 = -1.0;
    for (iterator di = begin(); di != end() && dui != dtu.end(); ++di)
    {
	if ( (*di - d0) > eps_dd )
	{
	    *(dui++) = *di;
	    d0 = *di;
	}
    }
    return dtu;
}

void DistanceTable::init()
{
    if (size() == 0)
    {
	NAtoms = 1;
	Nuniqd = 0;
	max_d = 1.0;
	return;
    }
    // calculate NAtoms
    double xNAtoms = 0.5 + sqrt(1 + 8.0*size())/2.0;
    NAtoms = int(xNAtoms);
    // sort and check values
    sort(begin(), end());
    if (at(0) <= 0)
    {
	cerr << "E: non-positive entry in DistanceTable, " <<
	    "d[0]=" << at(0) << endl;
	throw InvalidDistanceTable();
    }
    // calculate Nuniqd
    double eps_dd = sqrt(BGA::eps_badness);
    Nuniqd = 1;
    // here size() > 0
    for (iterator ilo = begin(), ihi = begin()+1; ihi != end(); ++ilo, ++ihi)
    {
	if ( (*ihi - *ilo) > eps_dd )
	    ++Nuniqd;
    }
    // max_d is very simple
    max_d = back();
}


////////////////////////////////////////////////////////////////////////
// Molecule definitions
////////////////////////////////////////////////////////////////////////

// static members
double Molecule::tol_dd  = numeric_limits<double>().max();
double Molecule::tol_nbad  = 0.05*0.05;
double Molecule::tol_r = 1.0e-8;
double Molecule::evolve_frac = 0.1;
bool Molecule::evolve_jump = true;
bool Molecule::evolve_relax = false;
bool Molecule::degenerate_relax = false;
vector<AtomFilter_t*> Molecule::atom_filters;
double Molecule::lookout_prob = 0.0;
Molecule::file_fmt_type Molecule::output_format = XYZ;

Molecule::Molecule()
{
    init();
}

Molecule::Molecule(const DistanceTable& dtab) : dTarget(dtab)
{
    init();
}

Molecule::Molecule(const DistanceTable& dtab,
	const int s, const double* px, const double* py, const double* pz
	) : dTarget(dtab)
{
    init();
    for (int i = 0; i < s; ++i)
    {
	Add(px[i], py[i], pz[i]);
    }
}

Molecule::Molecule(const DistanceTable& dtab,
	const vector<double>& vx, const vector<double>& vy,
	const vector<double>& vz
	) : dTarget(dtab)
{
    init();
    if (vx.size() != vy.size() || vx.size() != vz.size())
    {
	cerr << "E: invalid coordinate vectors" << endl;
	throw InvalidMolecule();
    }
    for (size_t i = 0; i != vx.size(); ++i)
    {
	Add(vx[i], vy[i], vz[i]);
    }
}

Molecule::Molecule(const Molecule& M) : dTarget(M.dTarget)
{
    init();
    *this  = M;
}

Molecule& Molecule::operator=(const Molecule& M)
{
    if (this == &M) return *this;
    // Clear() must be the first statement
    Clear();
    dTarget = M.dTarget;
    // duplicate source atoms
    atoms.resize(M.atoms.size(), NULL);
    for (AtomSequenceIndex seq(&M); !seq.finished(); seq.next())
    {
	atoms[seq.idx()] = new Atom_t(seq.ref());
    }
    pmx_used_distances = M.pmx_used_distances;
    pmx_partial_costs = M.pmx_partial_costs;
    free_pmx_slots = M.free_pmx_slots;
    // finished duplication
    max_natoms = M.max_natoms;
    badness = M.badness;
    // IO helpers
    output_format = M.output_format;
    opened_file = M.opened_file;
    trace = M.trace;
    return *this;
}

void Molecule::init()
{
    badness = 0;
    setMaxNAtoms(dTarget.NAtoms);
}

Molecule::~Molecule()
{
    // we must call Clear() to delete Atom_t objects
    Clear();
}


//////////////////////////////////////////////////////////////////////////
//// Molecule badness/fitness evaluation
//////////////////////////////////////////////////////////////////////////

void Molecule::Recalculate()
{
    if (NAtoms() > maxNAtoms())
    {
	cerr << "E: molecule too large in Recalculate()" << endl;
	throw InvalidMolecule();
    }
    // reset molecule
    badness = 0;
    // reset all atoms
    typedef vector<Atom_t*>::iterator VPAit;
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	seq.ptr()->ResetBadness();
    }
    // create array of all badness contributions with atom index
    typedef pair<double,int> BadnessWithIndex;
    vector<BadnessWithIndex> bwi;
    bwi.reserve(2*NDist());
    for (AtomSequenceIndex seq0(this); !seq0.finished(); seq0.next())
    {
	AtomSequenceIndex seq1 = seq0;
	for (seq1.next(); !seq1.finished(); seq1.next())
	{
	    int idx0 = seq0.ptr()->pmxidx;
	    int idx1 = seq1.ptr()->pmxidx;
	    double badnesshalf = pmx_partial_costs(idx0,idx1) / 2.0;
	    bwi.push_back(BadnessWithIndex(badnesshalf, seq0.idx()));
	    bwi.push_back(BadnessWithIndex(badnesshalf, seq1.idx()));
	}
    }
    // order pair iterators by corresponding badness for accurate summation
    sort(bwi.begin(), bwi.end());
    // sum over sorted errors
    for (vector<BadnessWithIndex>::iterator pbwi = bwi.begin();
	    pbwi != bwi.end(); ++pbwi)
    {
	atoms[pbwi->second]->IncBadness(pbwi->first);
	badness += pbwi->first;
    }
}

//void Molecule::ReassignPairs()

double Molecule::Badness() const
{
    return badness;
}

double Molecule::NormBadness() const
{
    return NDist() == 0 ? 0.0 : Badness()/NDist();
}

bool comp_Atom_Badness(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.Badness() < rhs.Badness();
}

bool comp_pAtom_Badness(const Atom_t* lhs, const Atom_t* rhs)
{
    return lhs->Badness() < rhs->Badness();
}


bool comp_pAtom_FreeBadness(const Atom_t* lhs, const Atom_t* rhs)
{
    return lhs->FreeBadness() < rhs->FreeBadness();
}


//////////////////////////////////////////////////////////////////////////
// Molecule operators
//////////////////////////////////////////////////////////////////////////

bool operator==(const Molecule& m0, const Molecule& m1)
{
    if (&m0 == &m1)
	return true;
    if (m0.maxNAtoms() != m1.maxNAtoms())
	return false;
    AtomSequence seq0(&m0), seq1(&m1);
    for (; !seq0.finished() && !seq1.finished() && *seq0.ptr() == *seq1.ptr();
	    seq0.next(), seq1.next() )
    { }
    return seq0.finished() && seq1.finished();
}

void Molecule::setMaxNAtoms(int sz)
{
    if (sz > dTarget.NAtoms && tol_dd > 0.0)
    {
	cerr << "E: not enough distances for maxNAtoms = " << sz << '.' <<
	    "  Did you forget tol_dd = 0?" << endl;
	throw InvalidMolecule();
    }
    else if (sz < 1)
    {
	cerr << "E: invalid value of maxNAtoms = " << sz << endl;
	throw InvalidMolecule();
    }
    else if (sz < NAtoms())
    {
	cerr << "E: molecule too large in setMaxNAtoms()" << endl;
	throw InvalidMolecule();
    }
    max_natoms = sz;
}

void Molecule::Shift(double dx, double dy, double dz)
{
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	Atom_t* pa = seq.ptr();
	pa->r[0] += dx;
	pa->r[1] += dy;
	pa->r[2] += dz;
    }
}

void Molecule::Center()
{
    double avg_rx = 0.0, avg_ry = 0.0, avg_rz = 0.0;
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	Atom_t* pa = seq.ptr();
	avg_rx += pa->r[0];
	avg_ry += pa->r[1];
	avg_rz += pa->r[2];
    }
    avg_rx /= NAtoms();
    avg_ry /= NAtoms();
    avg_rz /= NAtoms();
    Shift(-avg_rx, -avg_ry, -avg_rz);
}

void Molecule::Pop(const int aidx)
{
    if (aidx < 0 || aidx >= NAtoms())
    {
	throw range_error("in Molecule::Pop(const int aidx)");
    }
    // Pop should never get called on fixed atom
    assert(!atoms[aidx]->fixed);
    Atom_t* pa = atoms[aidx];
    removeAtomPairs(pa);
    free_pmx_slots.insert(pa->pmxidx);
    delete pa;
    atoms.erase(atoms.begin() + aidx);
}

void Molecule::Pop(const list<int>& cidx)
{
    // create a sorted set of indices of atoms to be popped
    set<int> popped(cidx.begin(), cidx.end());
    set<int>::reverse_iterator rii;
    for (rii = popped.rbegin(); rii != popped.rend(); ++rii)  Pop(*rii);
}

void Molecule::Clear()
{
    // return used distances
    for (AtomSequence seq0(this); !seq0.finished(); seq0.next())
    {
	AtomSequence seq1 = seq0;
	for (seq1.next(); !seq1.finished(); seq1.next())
	{
	    int i0 = seq0.ptr()->pmxidx;
	    int i1 = seq1.ptr()->pmxidx;
	    double& udst = pmx_used_distances(i0, i1);
	    if (udst > 0.0)	dTarget.push_back(udst);
	    udst = 0.0;
	}
    }
    sort(dTarget.begin(), dTarget.end());
    // remove all atoms
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	delete seq.ptr();
    }
    atoms.clear();
    free_pmx_slots.clear();
    badness = 0.0;
}

void Molecule::Add(const Molecule& M)
{
    for (AtomSequence seq(&M); !seq.finished(); seq.next())
    {
	Add(seq.ref());
    }
}

void Molecule::Add(double rx0, double ry0, double rz0)
{
    Add(Atom_t(rx0, ry0, rz0));
}

void Molecule::Add(const Atom_t& atom)
{
    if (NAtoms() == maxNAtoms())
    {
	cerr << "E: molecule too large in Add()" << endl;
	throw InvalidMolecule();
    }
    // create new atom
    Atom_t* pnew_atom;
    pnew_atom = new Atom_t(atom);
    pnew_atom->ResetBadness();
    pnew_atom->pmxidx = getPairMatrixIndex();
    // create new pairs while summing up the costs
    addNewAtomPairs(pnew_atom);
    atoms.push_back(pnew_atom);
}

void Molecule::Fix(const int cidx)
{
    if (cidx < 0 || cidx >= NAtoms())
    {
	ostringstream emsg;
	emsg << "Molecule::Fix(const int) invalid index " << cidx;
	throw range_error(emsg.str());
    }
    atoms[cidx]->fixed = true;
}

inline bool pAtom_is_fixed(const Atom_t* pa) { return pa->fixed; }

int Molecule::NFixed() const
{
    return count_if(atoms.begin(), atoms.end(), pAtom_is_fixed);
}

valarray<int> Molecule::good_neighbors_count(const vector<Atom_t>& vta)
{
    valarray<int> cnt(0, vta.size());
    // high limit for badness of a pair with good neighbor
    double hi_pbad = tol_nbad / 10.0;
    typedef vector<Atom_t>::const_iterator VAcit;
    VAcit a1i, a2i;
    int* pcnt1; int* pcnt2;
    for (  a1i = vta.begin(), pcnt1 = &(cnt[0]);
	    a1i != vta.end(); ++a1i, ++pcnt1 )
    {
	for (   a2i = a1i+1, pcnt2 = pcnt1+1; a2i != vta.end();
		++a2i, ++pcnt2 )
	{
	    double d = dist(*a1i, *a2i);
	    double dd = *(dTarget.find_nearest(d)) - d;
	    if (penalty(dd) < hi_pbad)
	    {
		++(*pcnt1);
		++(*pcnt2);
	    }
	}
    }
    return cnt;
}

void Molecule::filter_good_atoms(vector<Atom_t>& vta,
	double evolve_range, double hi_abad)
{
    if (NAtoms() == maxNAtoms())
    {
	cerr << "E: Molecule too large in filter_good_atoms()" << endl;
	throw InvalidMolecule();
    }
    typedef vector<Atom_t>::iterator VAit;
    VAit gai;
    // first check if atoms pass through atom_filters
    gai = vta.begin();
    if ( !atom_filters.empty() )
    {
	for (VAit tai = vta.begin(); tai != vta.end(); ++tai)
	{
	    if (check_atom_filters(&(*tai)))
	    {
		*(gai++) = *tai;
	    }
	}
	vta.erase(gai, vta.end());
    }
    // obtain badness of every test atom, do lazy evaluation when badness
    // is higher than (minimum + evolve_range)
    AtomCost* atomcost = getAtomCostCalculator();
    atomcost->setCutoff(hi_abad);
    atomcost->setCutoffRange(evolve_range);
    for (VAit tai = vta.begin(); tai != vta.end(); ++tai)
    {
	tai->IncBadness(atomcost->eval(*tai));
    }
    // atom cost cutoff is available here,
    // let us keep only good atoms
    gai = vta.begin();
    for (VAit tai = vta.begin(); tai != vta.end(); ++tai)
    {
	if (tai->Badness() > atomcost->cutoff())    continue;
	*(gai++) = *tai;
    }
    vta.erase(gai, vta.end());
}

bool Molecule::check_atom_filters(Atom_t* pa)
{
    typedef vector<AtomFilter_t*>::iterator VPAFit;
    bool isgood = true;
    for (   VPAFit pafi = atom_filters.begin();
	    isgood && pafi != atom_filters.end(); ++pafi )
    {
	isgood = (*pafi)->Check(pa, this);
    }
    return isgood;
}

struct rxa_par
{
    typedef vector<Atom_t*> VPA;
    typedef valarray<double> VAD;
    VPA* atoms;
    VAD* ad0;
    VAD* wt;
};

int rxa_fdf(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
{
    static Atom_t ta(0.0, 0.0, 0.0);
    copy(x->data, x->data + x->size, ta.r.data());
    AtomCost* atomcost = static_cast<AtomCost*>(params);
    atomcost->eval(ta);
    for (size_t m = 0; m != atomcost->lsqComponentsSize(); ++m)
    {
	gsl_vector_set(f, m, atomcost->lsqComponents()[m]);
	for (size_t n = 0; n != atomcost->lsqParametersSize(); ++n)
	{
	    gsl_matrix_set(J, m, n, atomcost->lsqJacobianGet(m, n));
	}
    }
    return GSL_SUCCESS;
}

int rxa_f(const gsl_vector* x, void* params, gsl_vector* f)
{
    static Atom_t ta(0.0, 0.0, 0.0);
    copy(x->data, x->data + x->size, ta.r.data());
    AtomCost* atomcost = static_cast<AtomCost*>(params);
    atomcost->eval(ta);
    for (size_t m = 0; m != atomcost->lsqComponentsSize(); ++m)
    {
	gsl_vector_set(f, m, atomcost->lsqComponents()[m]);
    }
    return GSL_SUCCESS;
}

int rxa_df(const gsl_vector* x, void* params, gsl_matrix* J)
{
    static Atom_t ta(0.0, 0.0, 0.0);
    copy(x->data, x->data + x->size, ta.r.data());
    AtomCost* atomcost = static_cast<AtomCost*>(params);
    atomcost->eval(ta);
    for (size_t m = 0; m != atomcost->lsqComponentsSize(); ++m)
    {
	for (size_t n = 0; n != atomcost->lsqParametersSize(); ++n)
	{
	    gsl_matrix_set(J, m, n, atomcost->lsqJacobianGet(m, n));
	}
    }
    return GSL_SUCCESS;
}

void Molecule::RelaxAtom(vector<Atom_t*>::iterator pai)
{
    RelaxAtom(pai - atoms.begin());
}

void Molecule::RelaxAtom(const int cidx)
{
    if (cidx < 0 || cidx >= NAtoms())
    {
	throw range_error("in Molecule::RelaxAtom(const int)");
    }
    // RelaxAtom should never get called on a fixed atom
    assert(!atoms[cidx]->fixed);
    Atom_t ta = getAtom(cidx);
    Pop(cidx);
    RelaxExternalAtom(ta);
    Add(ta);
}

void Molecule::RelaxExternalAtom(Atom_t& ta)
{
    const int maximum_relaxations = 20;
    const int maximum_iterations = 500;
    // pj: this seems to be crashing when NAtoms() < 3
    if (NAtoms() < 3)
	return;
    // do relaxation on a copy of ta
    Atom_t rta(ta);
    // loop while badness is improved
    double lo_abad = DOUBLE_MAX;
    AtomCost* atomcost = getAtomCostCalculator();
    for (int nrelax = 0; nrelax < maximum_relaxations; ++nrelax)
    {
	double tbad = atomcost->eval(rta);
	// get out if tbad did not improve
	if (!eps_lt(tbad, lo_abad))	break;
	lo_abad = tbad;
	ta = rta;
	// get out if lo_abad is very low
	if (lo_abad < BGA::eps_badness)	break;
	// carry out relaxation otherwise
	// define function to be minimized
	gsl_multifit_function_fdf f;
	f.f = &rxa_f;
	f.df = &rxa_df;
	f.fdf = &rxa_fdf;
	f.n = atomcost->lsqComponentsSize();
	f.p = atomcost->lsqParametersSize();
	f.params = static_cast<void*>(atomcost);
	// bind rta coordinates to vector x
	gsl_vector_view x = gsl_vector_view_array(rta.r.data(), 3);
	// allocate solver
	gsl_multifit_fdfsolver* lms = gsl_multifit_fdfsolver_alloc(
		gsl_multifit_fdfsolver_lmsder, f.n, f.p);
	gsl_multifit_fdfsolver_set(lms, &f, &x.vector);
	// allocate vector for gradient test
	gsl_vector* G = gsl_vector_alloc(3);
	// minimize atom badness
	int iter = 0, status;
	do
	{
	    ++iter;
	    status = gsl_multifit_fdfsolver_iterate(lms);
	    if (status)
	    {
		if ( status != GSL_ETOLF && status != GSL_ETOLX &&
			status != GSL_CONTINUE )
		{
		    cerr << "LM solver status = " << status << " " <<
			gsl_strerror(status) << endl;
		}
		break;
	    }
	    // scale f_i with atom fitness
	    // pj: test gradient or do a simplex search?
	    gsl_multifit_gradient(lms->J, lms->f, G);
	    status = gsl_multifit_test_gradient(G, BGA::eps_badness/tol_r);
//	    status = gsl_multifit_test_delta(lms->dx, lms->x, tol_r, tol_r);
	}
	while (status == GSL_CONTINUE && iter < maximum_iterations);
	for (int i = 0; i < 3; ++i)
	    rta.r[i] = gsl_vector_get(lms->x, i);
	gsl_vector_free(G);
	gsl_multifit_fdfsolver_free(lms);
    }
}

AtomCost* Molecule::getAtomCostCalculator()
{
    static AtomCost the_acc(this);
    the_acc.resetFor(this);
    return &the_acc;
}

void Molecule::addNewAtomPairs(Atom_t* pa)
{
    // calculate the cost
    AtomCost* atomcost = getAtomCostCalculator();
    atomcost->eval(pa);
    // store partial costs
    const vector<double>& ptcs = atomcost->partialCosts();
    assert(atoms.size() == ptcs.size());
    for (AtomSequenceIndex seq(this); !seq.finished(); seq.next())
    {
	double pairbadness = ptcs[seq.idx()];
	int idx0 = pa->pmxidx;
	int idx1 = seq.ptr()->pmxidx;
	pmx_partial_costs(idx0,idx1) = pairbadness;
	double badnesshalf = pairbadness / 2.0;
	seq.ptr()->IncBadness(badnesshalf);
	pa->IncBadness(badnesshalf);
    }
    badness += atomcost->total();
    if (badness < BGA::eps_badness)	badness = 0.0;    
    // remember used distances in pmx_used_distances
    const vector<int>& didcs = atomcost->usedTargetDistanceIndices();
    const vector<int>& aidcs = atomcost->usedTargetAtomIndices();
    assert(didcs.size() == aidcs.size());
    vector<int>::const_iterator dii = didcs.begin();
    vector<int>::const_iterator aii = aidcs.begin();
    for (; dii != didcs.end(); ++dii, ++aii)
    {
	int idx0 = pa->pmxidx;
	int idx1 = atoms[*aii]->pmxidx;
	pmx_used_distances(idx0,idx1) = dTarget[*dii];
    }
    // remove used distances from target table
    set<int> uidcs(didcs.begin(), didcs.end());
    set<int>::reverse_iterator ii;
    for (ii = uidcs.rbegin(); ii != uidcs.rend(); ++ii)
    {
	dTarget.erase(dTarget.begin() + *ii);
    }
}

void Molecule::removeAtomPairs(Atom_t* pa)
{
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	if (pa == seq.ptr())	continue;
	// return any used distances
	int idx0 = pa->pmxidx;
	int idx1 = seq.ptr()->pmxidx;
	double& udst = pmx_used_distances(idx0, idx1);
	if (udst > 0.0)	    dTarget.return_back(udst);
	udst = 0.0;
	// remove pair costs
	double pairbadness = pmx_partial_costs(idx0,idx1);
	double badnesshalf = pairbadness/2.0;
	pa->DecBadness(badnesshalf);
	seq.ptr()->DecBadness(badnesshalf);
	badness -= pairbadness;
    }
    if (badness < BGA::eps_badness)	badness = 0.0;    
}

int Molecule::push_good_distances(
	vector<Atom_t>& vta, double* afit, int ntrials
	)
{
    if (!ntrials)   return 0;
    // add new atom in direction defined by 2 atoms
    if (NAtoms() == maxNAtoms())
    {
	cerr << "E: molecule too large for finding a new position" << endl;
	throw InvalidMolecule();
    }
    else if (NAtoms() < 1)
    {
	cerr << "E: empty molecule, no way to push_good_distances()" << endl;
	throw InvalidMolecule();
    }
    const double eps_d = 10.0*sqrt(numeric_limits<double>().epsilon());
    int push_count = 0;
    for (int nt = 0; nt < ntrials; ++nt)
    {
        size_t basesize = min(NAtoms(),2);
	const PickType& aidx = randomWeighedPick(basesize, NAtoms(), afit);
	Atom_t& a0 = *atoms[aidx[0]];
        R3::Vector rdir(0.0, 0.0, 0.0);
	if (NAtoms() > 1)
	{
	    Atom_t& a1 = *atoms[aidx[1]];
            rdir = a1.r - a0.r;
	}
	// normalize rdir if defined
	bool lattice_rdir;
	double nm_rdir = R3::norm(rdir);
	if (nm_rdir > eps_d)
	{
	    rdir /= nm_rdir;
	    lattice_rdir = true;
	}
	// otherwise orient along the z-axis
	else
	{
	    rdir[0] = 0.0;
	    rdir[1] = 0.0;
	    rdir[2] = 1.0;
	    lattice_rdir = false;
	}
	// pick free distance
	int didx = randomInt(dTarget.size());
	double radius = dTarget[didx];
	// add front atom
        R3::Vector nr;
        nr = a0.r + rdir*radius;
	Atom_t ad1front(nr);
	ad1front.ttp = LINEAR;
	vta.push_back(ad1front);
	++push_count;
	// check opposite direction when it makes sense
	// this accounts for extra trial
	if (lattice_rdir)
	{
	    ++nt;
            nr = a0.r - rdir*radius;
	    Atom_t ad1back(nr);
	    ad1back.ttp = LINEAR;
	    vta.push_back(ad1back);
	    ++push_count;
	}
    }
    return push_count;
}

int Molecule::push_good_triangles(
	vector<Atom_t>& vta, double* afit, int ntrials
	)
{
    if (!ntrials)   return 0;
    // generate randomly oriented triangles
    if (NAtoms() == maxNAtoms())
    {
	cerr << "E: molecule too large for finding a new position" << endl;
	throw InvalidMolecule();
    }
    else if (NAtoms() < 2)
    {
	cerr << "E: molecule too small, triangulation not possible" << endl;
	throw InvalidMolecule();
    }
    const double eps_d = 10.0*sqrt(numeric_limits<double>().epsilon());
    int push_count = 0;
    for (int nt = 0; nt < ntrials; ++nt)
    {
	// pick 2 atoms for base and 3rd for plane orientation
	int nchoose = NAtoms() > 2 ? 3 : 2;
	const PickType& aidx = randomWeighedPick(nchoose, NAtoms(), afit);
	Atom_t& a0 = *atoms[aidx[0]];
	Atom_t& a1 = *atoms[aidx[1]];
	// pick 2 vertex distances
	bool with_repeat = (tol_dd <= 0.0);
        const PickType& didx = with_repeat ?
            randomPickWithRepeat(2, dTarget.size()) :
            randomPickFew(2, dTarget.size());
	double r13 = dTarget[didx[0]];
	double r23 = dTarget[didx[1]];
	double r12 = dist(a0, a1);
	// is triangle base reasonably large?
	if (r12 < eps_d)
	    continue;
	// get and store both possible values of xlong
	double xl0 = (r13*r13 + r12*r12 - r23*r23) / (2.0*r12);
	double xlong[2] = { xl0, r12-xl0 };
	// get and store both possible values of xperp
	double xp2 = r13*r13 - xlong[0]*xlong[0];
	double xp = sqrt(fabs(xp2));
	if (xp < eps_d)
	    xp = 0.0;
	else if (xp2 < 0.0)
	    continue;
	double xperp[2] = { -xp, xp };
	// find direction along triangle base:
        R3::Vector longdir;
        longdir = (a1.r - a0.r)/r12;
	// generate direction perpendicular to longdir
        R3::Vector perpdir(0.0, 0.0, 0.0);
	if (nchoose > 2)
	{
	    Atom_t& a2 = *atoms[aidx[2]];
            perpdir = a2.r - a0.r;
	    perpdir -= longdir*R3::dot(longdir, perpdir);
	}
	// normalize perpdir if defined
	bool lattice_plane;
	double nm_perpdir = R3::norm(perpdir);
	if (nm_perpdir > eps_d)
	{
	    perpdir /= nm_perpdir;
	    lattice_plane = true;
	}
	// otherwise generate direction in cartesian plane
	// perpendicular to the smallest component of longdir
	else
	{
            R3::Vector uv(0.0, 0.0, 0.0);
            R3::Vector ald = fabs(longdir);
	    int ijk = minIndex(ald);
	    uv[ijk] = 1.0;
	    perpdir = R3::cross(longdir, uv);
	    perpdir /= R3::norm(perpdir);
	    lattice_plane = false;
	}
	// allocate vallarays for positions of a0 and vertex P
        R3::Vector Pa0 = a0.r;
        R3::Vector P;
	// if vertex search has already failed above, nt would increase by 1
	// here we want nt to count number of added vertices
	--nt;
	// loops over all 4 vertices in case of lattice_plane
	for (double* pxl = xlong; pxl != xlong+2; ++pxl)
	{
	    for (double* pxp = xperp; pxp != xperp+2; ++pxp)
	    {
		++nt;
		P = Pa0 + (*pxl)*longdir + (*pxp)*perpdir;
		Atom_t ad2(P);
		ad2.ttp = PLANAR;
		vta.push_back(ad2);
		++push_count;
		if (!lattice_plane)
		    break;
	    }
	    if (!lattice_plane)
		break;
	}
    }
    return push_count;
}

int Molecule::push_good_pyramids(
	vector<Atom_t>& vta, double* afit, int ntrials
	)
{
    if (!ntrials)   return 0;
    if (NAtoms() == maxNAtoms())
    {
	cerr << "E: molecule too large for finding a new position" << endl;
	throw InvalidMolecule();
    }
    else if (NAtoms() < 3)
    {
	cerr << "E: molecule too small, cannot construct pyramid" << endl;
	throw InvalidMolecule();
    }
    const double eps_d = 10.0*sqrt(numeric_limits<double>().epsilon());
    int push_count = 0;
    for (int nt = 0; nt < ntrials;)
    {
	// pick 3 base atoms
        const PickType& aidx = randomWeighedPick(3, NAtoms(), afit);
	Atom_t& a0 = *atoms[aidx[0]];
	Atom_t& a1 = *atoms[aidx[1]];
	Atom_t& a2 = *atoms[aidx[2]];
	// pick 3 vertex distances
	bool with_repeat = (tol_dd <= 0.0);
	PickType didx = with_repeat ?
            randomPickWithRepeat(3, dTarget.size()) :
            randomPickFew(3, dTarget.size());
	// loop over all permutations of selected distances
	sort(didx.begin(), didx.end());
	do
	{
	    ++nt;
	    double r14 = dTarget[ didx[0] ];
	    double r24 = dTarget[ didx[1] ];
	    double r34 = dTarget[ didx[2] ];
	    // uvi is a unit vector in a0a1 direction
            R3::Vector uvi;
            uvi = a1.r - a0.r;
	    double r12 = R3::norm(uvi);
	    if (r12 < eps_d)
		continue;
	    uvi /= r12;
	    // v13 is a0a2 vector
            R3::Vector v13;
            v13 = a2.r - a0.r;
	    // uvj lies in a0a1a2 plane and is perpendicular to uvi
            R3::Vector uvj = v13;
	    uvj -= uvi*R3::dot(uvi, uvj);
	    double nm_uvj = R3::norm(uvj);
	    if (nm_uvj < eps_d)
		continue;
	    uvj /= nm_uvj;
	    // uvk is a unit vector perpendicular to a0a1a2 plane
            R3::Vector uvk = R3::cross(uvi, uvj);
	    double xP1 = -0.5/(r12)*(r12*r12 + r14*r14 - r24*r24);
	    // Pn are coordinates in pyramid coordinate system
            R3::Vector P1;
	    P1[0] = xP1; P1[1] = P1[2] = 0.0;
	    // vT is translation from pyramid to normal system
            R3::Vector vT;
	    vT[0] = a0.r[0] - xP1*uvi[0];
	    vT[1] = a0.r[1] - xP1*uvi[1];
	    vT[2] = a0.r[2] - xP1*uvi[2];
	    // obtain coordinates of P3
            R3::Vector P3 = P1;
	    P3[0] += R3::dot(uvi, v13);
	    P3[1] += R3::dot(uvj, v13);
	    P3[2] = 0.0;
	    double xP3 = P3[0];
	    double yP3 = P3[1];
	    // find pyramid vertices
            R3::Vector P4(0.0, 3);
	    double h2 = r14*r14 - xP1*xP1;
	    // does P4 belong to a0a1 line?
	    if (fabs(h2) < eps_d)
	    {
		// is vertex on a0a1
		if (fabs(R3::norm(P3) - r14) > eps_d)
		    continue;
		P4 = vT;
		Atom_t ad3(P4[0], P4[1], P4[2]);
		ad3.ttp = SPATIAL;
		vta.push_back(ad3);
		++push_count;
		continue;
	    }
	    else if (h2 < 0)
	    {
		continue;
	    }
	    double yP4 = 0.5/(yP3)*(h2 + xP3*xP3 + yP3*yP3 - r34*r34);
	    double z2P4 = h2 - yP4*yP4;
	    // does P4 belong to a0a1a2 plane?
	    if (fabs(z2P4) < eps_d)
	    {
		P4 = yP4*uvj + vT;
		Atom_t ad3(P4[0], P4[1], P4[2]);
		ad3.ttp = SPATIAL;
		vta.push_back(ad3);
		++push_count;
		continue;
	    }
	    else if (z2P4 < 0)
	    {
		continue;
	    }
	    // here we can construct 2 pyramids
	    double zP4 = sqrt(z2P4);
	    // top one
	    P4 = yP4*uvj + zP4*uvk + vT;
	    Atom_t ad3top(P4[0], P4[1], P4[2]);
	    ad3top.ttp = SPATIAL;
	    vta.push_back(ad3top);
	    ++push_count;
	    // and bottom one, which makes an extra trial
	    ++nt;
	    P4 = yP4*uvj - zP4*uvk + vT;
	    Atom_t ad3bottom(P4[0], P4[1], P4[2]);
	    ad3bottom.ttp = SPATIAL;
	    vta.push_back(ad3bottom);
	    push_count++;
	} while ( next_permutation(didx.begin(), didx.end()) );
    }
    return push_count;
}

int Molecule::push_second_atoms(vector<Atom_t>& vta, int ntrials)
{
    if (NAtoms() != 1)
    {
	cerr << "E: push_second_atoms() must be called with 1-atom molecule"
	    << endl;
	throw InvalidMolecule();
    }
    // second atoms will be pushed along z-direction from a0
    Atom_t& a0 = *atoms[0];
    // new position
    R3::Vector nr = a0.r;
    int push_count = 0;
    if (ntrials > 2*dTarget.Nuniqd)
    {
	// we can push all the unique distances in both directions
	typedef vector<double>::iterator VDit;
	vector<double> dtu = dTarget.unique();
	for (VDit dui = dtu.begin(); dui != dtu.end(); ++dui)
	{
	    // top atom
	    nr[2] = a0.r[2] + *dui;
	    Atom_t adtop(nr); adtop.ttp = LINEAR;
	    vta.push_back(adtop);
	    // bottom atom
	    nr[2] = a0.r[2] - *dui;
	    Atom_t adbot(nr); adbot.ttp = LINEAR;
	    vta.push_back(adbot);
	    push_count += 2;
	}
    }
    else
    {
	// distances will be picked randomly
	for (push_count = 0; push_count < ntrials; ++push_count)
	{
	    int didx = randomInt(dTarget.size());
	    double dz = dTarget[didx];
	    // randomize direction
            dz *= plusminus();
	    nr[2] = a0.r[2] + dz;
	    Atom_t ad(nr); ad.ttp = LINEAR;
	    vta.push_back(ad);
	}
    }
    return push_count;
}

int Molecule::push_third_atoms(vector<Atom_t>& vta, int ntrials)
{
    if (NAtoms() != 2)
    {
	cerr << "E: push_third_atoms() must be called with 2-atom molecule"
	    << endl;
	throw InvalidMolecule();
    }
    // build list of distances from the base atoms
    list<double> d0, d1;
    if (ntrials > 2 * dTarget.Nuniqd * dTarget.Nuniqd)
    {
	// we can push all the unique triangles
	typedef vector<double>::iterator VDit;
	vector<double> dtu = dTarget.unique();
	for (VDit ui0 = dtu.begin(); ui0 != dtu.end(); ++ui0)
	{
	    for (VDit ui1 = dtu.begin(); ui1 != dtu.end(); ++ui1)
	    {
		d0.push_back(*ui0);
		d1.push_back(*ui1);
	    }
	}
    }
    else
    {
	// distances will be picked randomly
	bool with_repeat = (tol_dd <= 0.0);
	for (int i = 0; i < ntrials; ++i)
	{
            const PickType& didx = with_repeat ?
                randomPickWithRepeat(2, dTarget.size()) :
                randomPickFew(2, dTarget.size());
	    d0.push_back(dTarget[didx[0]]);
	    d1.push_back(dTarget[didx[1]]);
	}
    }
    // define some constants and base atoms
    const double eps_d = 10.0*sqrt(numeric_limits<double>().epsilon());
    int push_count = 0;
    Atom_t a0 = *atoms[0];
    Atom_t a1 = *atoms[1];
    double r01 = dist(a0, a1);
    // longitudinal direction along triangle base
    R3::Vector longdir;
    longdir = (a1.r - a0.r)/r01;
    // perpendicular direction to longdir is a cross product of x with longdir
    R3::Vector ux(0.0, 0.0, 0.0);
    ux[0] = 1.0;
    R3::Vector perpdir = R3::cross(ux, longdir);
    double nm_perpdir = R3::norm(perpdir);
    if (nm_perpdir == 0.0)
    {
	// longdir is ux, let us set perpdir = uy
	perpdir = 0.0;
	perpdir[1] = 1.0;
    }
    else
    {
	perpdir /= nm_perpdir;
    }
    R3::Vector Pa0 = a0.r;
    // now loop over list of distances
    typedef list<double>::iterator LDit;
    for (   LDit d0i = d0.begin(), d1i = d1.begin();
	    d0i != d0.end() && d1i != d1.end(); ++d0i, ++d1i )
    {
	double& r02 = *d0i;
	double& r12 = *d1i;
	// is triangle base reasonably large?
	if (r01 < eps_d)
	    continue;
	// calculate xlong
	double xlong = (r02*r02 + r01*r01 - r12*r12) / (2.0*r01);
	// calculate xperp
	double xp2 = r02*r02 - xlong*xlong;
	double xperp = sqrt(fabs(xp2));
	if (xperp < eps_d)
	    xperp = 0.0;
	else if (xp2 < 0.0)
	    continue;
	else if (randomInt(2) == 0)
	    xperp = -xperp;
	// add atom
        R3::Vector Pn;
	Pn = Pa0 + xlong*longdir + xperp*perpdir;
	Atom_t ad(Pn); ad.ttp = PLANAR;
	vta.push_back(ad);
	++push_count;
    }
    return push_count;
}

void Molecule::Evolve(const int* est_triang)
{
    const int& nlinear = est_triang[LINEAR];
    const int& nplanar = est_triang[PLANAR];
    const int& nspatial = est_triang[SPATIAL];
    if (NAtoms() == maxNAtoms())
    {
	cerr << "E: full-sized molecule cannot Evolve()" << endl;
	throw InvalidMolecule();
    }
    // containter for test atoms
    vector<Atom_t> vta;
    // calculate array of atom fitnesses
    valarray<double> vafit(NAtoms());
    double* pd = &vafit[0];
    // first fill the array with badness
    typedef vector<Atom_t*>::iterator VPAit;
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	*(pd++) = seq.ptr()->Badness();
    }
    // finally get the reciprocal value
    vafit = recipw0(vafit);
    double* afit = &vafit[0];
    bool lookout = NAtoms() && NAtoms() <= 2 &&
	randomFloat() < lookout_prob;
    const int lookout_trials = 1500;
    // evolution is trivial for empty or 1-atom molecule
    switch (NAtoms())
    {
	case 0:
	    Add(0.0, 0.0, 0.0);
	    return;
	case 1:
	    if (lookout)
	    {
		push_second_atoms(vta, lookout_trials);
		break;
	    }
	case 2:
	    if (lookout)
	    {
		push_third_atoms(vta, lookout_trials);
		break;
	    }
	default:
	    push_good_distances(vta, afit, nlinear);
	    push_good_triangles(vta, afit, nplanar);
	    push_good_pyramids(vta, afit, nspatial);
    }
    // set badness range from desired badness
    double evolve_range = NAtoms()*tol_nbad*evolve_frac;
    double hi_abad = DOUBLE_MAX;
    // try to add as many atoms as possible
    typedef vector<Atom_t>::iterator VAit;
    while (true)
    {
	filter_good_atoms(vta, evolve_range, hi_abad);
	if (vta.size() == 0)
	    break;
	// calculate fitness of test atoms
	valarray<double> vtafit(vta.size());
	if (lookout)
	{
	    // calculate fitness as the number of good neighbors
	    valarray<int> cnt = good_neighbors_count(vta);
	    int max_cnt = cnt.max();
	    double* pfit = &vtafit[0];
	    int* pcnt = &cnt[0];
	    for ( ; pfit != &vtafit[vtafit.size()]; ++pfit, ++pcnt)
	    {
		*pfit = (*pcnt < max_cnt/2) ? 0.0 : *pcnt;
	    }
	}
	else
	{
	    // calculate fitness as reciprocal value of badness
	    // fill the vtafit array with badness
	    double* pfit = &vtafit[0];
	    for (VAit ai = vta.begin(); ai != vta.end(); ++ai, ++pfit)
	    {
		*pfit = ai->Badness();
	    }
	    // then get the reciprocal value
	    vtafit = recipw0(vtafit);
	}
	// vtafit is ready here
	int idx = randomWeighedInt(vtafit.size(), &vtafit[0]);
	Add(vta[idx]);
	hi_abad = vta[idx].Badness() + evolve_range;
	vta.erase(vta.begin()+idx);
	if (evolve_relax)
	{
	    VPAit worst = max_element(atoms.begin(), atoms.end(),
		    comp_pAtom_FreeBadness);
	    if ( eps_gt((*worst)->Badness(), 0.0) && ! (*worst)->fixed )
	    {
		RelaxAtom(worst);
	    }
	}
	if (NAtoms() == maxNAtoms() || !evolve_jump)
	    break;
	for (VAit pai = vta.begin(); pai != vta.end(); ++pai)
	    pai->ResetBadness();
    }
}

void Molecule::Degenerate(int Npop)
{
    Npop = min(NAtoms(), Npop);
    if (Npop == 0)  return;
    // build array of atom badnesses
    double freebad[NAtoms()];
    int freeidx[NAtoms()];
    int Nfree = 0;
    for (int i = 0; i != NAtoms(); ++i)
    {
	Atom_t* pai = atoms[i];
	if ( pai->fixed )  continue;
	freebad[Nfree] = pai->Badness();
	freeidx[Nfree] = i;
	Nfree++;
    }
    if (Nfree == 0)  return;
    Npop = min(Nfree, Npop);
    // build list of indices of atoms to pop
    const PickType& idxidx = randomWeighedPick(Npop, Nfree, freebad);
    list<int> ipop;
    PickType::const_iterator ii;
    for (ii = idxidx.begin(); ii != idxidx.end(); ++ii)
    {
	ipop.push_back(freeidx[*ii]);
    }
    Pop(ipop);
    typedef vector<Atom_t*>::iterator VPAit;
    if (degenerate_relax && NAtoms() > 1)
    {
	VPAit worst = max_element(atoms.begin(), atoms.end(),
		comp_pAtom_FreeBadness);
	if ( eps_gt((*worst)->Badness(), 0.0) && ! (*worst)->fixed )
	{
	    RelaxAtom(worst);
	}
    }
}

int Molecule::getPairMatrixIndex()
{
    int idx;
    if (!free_pmx_slots.empty())
    {
	set<int>::iterator firstfree = free_pmx_slots.begin();
	idx = *firstfree;
	free_pmx_slots.erase(firstfree);
	return idx;
    }
    // no free slots here
    idx = atoms.size();
    // check if we need to resize pmx_used_distances
    if (size_t(idx) >= pmx_used_distances.rows())
    {
	size_t sz0 = idx + 1;
	size_t sz1 = min(2*pmx_used_distances.rows(), size_t(max_natoms));
	size_t sz = max(sz0, sz1);
	assert(sz <= size_t(max_natoms));
	pmx_used_distances.resize(sz, sz, 0.0);
	pmx_partial_costs.resize(sz, sz, 0.0);
    }
    return idx;
}

////////////////////////////////////////////////////////////////////////
// Molecule IO functions
////////////////////////////////////////////////////////////////////////

Molecule::ParseHeader::ParseHeader(const string& s) : header(s)
{
    // parse format
    string fmt;
    // initialize members:
    state =
	read_token("BGA molecule format", fmt) &&
	read_token("NAtoms", NAtoms);
    if (!state)
    {
	if (read_token("Number of particles", NAtoms))
	{
	    format = ATOMEYE;
	    fmt = "atomeye";
	}
	else
	    return;
    }
    else if (fmt == "xyz")
	format = XYZ;
    else if (fmt == "atomeye")
	format = ATOMEYE;
    else
    {
	state = false;
	return;
    }
}

template<typename T> bool Molecule::ParseHeader::read_token(
	const char* token, T& value
	)
{
    const char* fieldsep = ":= ";
    int ltoken = strlen(token);
    string::size_type sp;
    const string::size_type npos = string::npos;
    if (
	    npos == (sp = header.find(token)) ||
	    npos == (sp = header.find_first_not_of(fieldsep, sp+ltoken))
       )
    {
	return false;
    }
    istringstream istrs(header.substr(sp));
    bool result = (istrs >> value);
    return result;
}

istream& Molecule::ReadXYZ(istream& fid)
{
    // read values to integer vector vxyz
    string header;
    vector<double> vxyz;
    bool result = read_header(fid, header) && read_data(fid, vxyz);
    if (!result) return fid;
    int vxyz_NAtoms = vxyz.size()/3;
    // check how many numbers were read
    ParseHeader ph(header);
    if (ph && vxyz_NAtoms != ph.NAtoms)
    {
	cerr << "E: " << opened_file << ": expected " << ph.NAtoms <<
	    " atoms, read " << vxyz_NAtoms << endl;
	throw IOError();
    }
    if ( vxyz.size() % 3 )
    {
	cerr << "E: " << opened_file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    for (vector<double>::iterator ii = vxyz.begin(); ii < vxyz.end();)
    {
	Add(Atom_t(*ii++, *ii++, *ii++));
    }
    return fid;
}

bool Molecule::ReadXYZ(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    opened_file = file;
    bool result = ReadXYZ(fid);
    opened_file.clear();
    fid.close();
    return result;
}

bool Molecule::WriteFile(const char* file)
{
    // test if file is writeable
    ofstream fid(file, ios_base::out|ios_base::ate );
    if (!fid)
    {
	ostringstream oss;
	oss << "WriteFile(): unable to write to '" << file << "'";
	cerr << oss.str() << endl;
	throw IOError(oss.str());
    }
    fid.close();
    // write via temporary file
    int filelen = strlen(file);
    char writefile[filelen+6+1];
    strcpy(writefile, file);
    memset(writefile+filelen, 'X', 6);
    writefile[filelen+6] = '\0';
    mktempofstream(fid, writefile);
    bool result = (fid << *this);
    fid.close();
    rename(writefile, file);
    return result;
}

bool Molecule::WriteXYZ(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtXYZ();
    bool result = WriteFile(file);
    output_format = org_ofmt;
    return result;
}

bool Molecule::WriteAtomEye(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtAtomEye();
    bool result = WriteFile(file);
    output_format = org_ofmt;
    return result;
}

void Molecule::OutFmtXYZ()
{
    output_format = XYZ;
}

void Molecule::OutFmtAtomEye()
{
    output_format = ATOMEYE;
}

istream& operator>>(istream& fid, Molecule& M)
{
    string header;
    istream::pos_type p = fid.tellg();
    if( !read_header(fid, header) )
    {
	fid.setstate(ios_base::failbit);
	return fid;
    }
    fid.seekg(p);
    Molecule::ParseHeader ph(header);
    if (!ph)
    {
	fid.setstate(ios_base::failbit);
	return fid;
    }
    bool result = true;
    switch (ph.format)
    {
	case Molecule::XYZ:
	    result = M.ReadXYZ(fid);
	    break;
	case Molecule::ATOMEYE:
	    throw runtime_error("reading of atomeye files not implemented");
	    break;
    }
    if (!result)
    {
	fid.setstate(ios_base::failbit);
    }
    return fid;
}

ostream& operator<<(ostream& fid, Molecule& M)
{
    typedef vector<Atom_t*>::iterator VPAit;
    VPAit afirst = M.atoms.begin();
    VPAit alast = M.atoms.end();
    switch (M.output_format)
    {
	case Molecule::XYZ:
	    for (VPAit pai = afirst; pai != alast; ++pai)
	    {
		fid <<
		    (*pai)->r[0] << " " <<
		    (*pai)->r[1] << " " <<
		    (*pai)->r[2] << endl;
	    }
	    break;
	case Molecule::ATOMEYE:
	    double xyz_lo = 0.0;
	    double xyz_hi = 1.0;
	    double xyz_range = xyz_hi - xyz_lo;
	    if (M.NAtoms() > 0)
	    {
		const double scale = 1.01;
		list<double> xyz_extremes;
		xyz_extremes.push_back(-M.dTarget.max_d);
		xyz_extremes.push_back(+M.dTarget.max_d);
		for (VPAit pai = afirst; pai != alast; ++pai)
		{
                    R3::Vector rsc = (*pai)->r * scale;
                    xyz_extremes.insert(xyz_extremes.end(),
                            rsc.data(), rsc.data() + 3);
		}
		// make atomeye happy
		xyz_extremes.push_back(-1.75);
		xyz_extremes.push_back(+1.75);
		xyz_lo = *min_element(xyz_extremes.begin(), xyz_extremes.end());
		xyz_hi = *max_element(xyz_extremes.begin(), xyz_extremes.end());
		xyz_range = xyz_hi - xyz_lo;
	    }
	    double xyz_med = (xyz_hi + xyz_lo)/2.0;
	    fid << "# BGA molecule format = atomeye" << endl;
	    fid << "# NAtoms = " << M.NAtoms() << endl;
	    fid << "Number of particles = " << M.NAtoms() << endl;
	    fid << "A = 1.0 Angstrom (basic length-scale)" << endl;
	    fid << "H0(1,1) = " << xyz_range << " A" << endl;
	    fid << "H0(1,2) = 0 A" << endl;
	    fid << "H0(1,3) = 0 A" << endl;
	    fid << "H0(2,1) = 0 A" << endl;
	    fid << "H0(2,2) = " << xyz_range << " A" << endl;
	    fid << "H0(2,3) = 0 A" << endl;
	    fid << "H0(3,1) = 0 A" << endl;
	    fid << "H0(3,2) = 0 A" << endl;
	    fid << "H0(3,3) = " << xyz_range << " A" << endl;
	    fid << ".NO_VELOCITY." << endl;
	    // 4 entries: x, y, z, Uiso
	    fid << "entry_count = 4" << endl;
	    fid << "auxiliary[0] = abad [au]" << endl;
	    fid << endl;
	    // pj: now it only works for a single Carbon atom in the molecule
	    fid << "12.0111" << endl;
	    fid << "C" << endl;
	    for (VPAit pai = afirst; pai != alast; ++pai)
	    {
		fid <<
		    ((*pai)->r[0] - xyz_med) / xyz_range + 0.5 << " " <<
		    ((*pai)->r[1] - xyz_med) / xyz_range + 0.5 << " " <<
		    ((*pai)->r[2] - xyz_med) / xyz_range + 0.5 << " " <<
		    (*pai)->Badness() << endl;
	    }
	    break;
    }
    return fid;
}

void Molecule::PrintBadness() const
{
    if (!NAtoms())  return;
    cout << "ABadness() =";
    double mab = (*max_element(atoms.begin(), atoms.end(),
		    comp_pAtom_Badness)) -> Badness();
    bool marked = false;
    vector<Atom_t*>::const_iterator pai;
    for (pai = atoms.begin(); pai != atoms.end(); ++pai)
    {
	cout << ' ';
	if (!marked && (*pai)->Badness() == mab)
	{
	    cout << '+';
	    marked = true;
	}
	cout << (*pai)->Badness();
    }
    cout << endl;
}

void Molecule::PrintFitness()
{
    valarray<double> vafit(NAtoms());
    double* pd = &vafit[0];
    typedef vector<Atom_t*>::iterator VPAit;
    // first fill the array with badness
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai, ++pd)
	*pd = (*pai)->Badness();
    // then get the reciprocal value
    vafit = recipw0(vafit);
    cout << "AFitness() =";
    double mab = vafit.max();
    bool marked = false;
    for (pd = &vafit[0]; pd != &vafit[vafit.size()]; ++pd)
    {
	cout << ' ';
	if (!marked && *pd == mab)
	{
	    cout << '+';
	    marked = true;
	}
	cout << *pd;
    }
    cout << endl;
}

// End of file
