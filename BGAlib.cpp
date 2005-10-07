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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "BGAutils.hpp"
#include "BGAlib.hpp"

// random number generator
gsl_rng * BGA::rng = gsl_rng_alloc(gsl_rng_default);
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


////////////////////////////////////////////////////////////////////////
// Atom_t definitions
////////////////////////////////////////////////////////////////////////

Atom_t::Atom_t(double rx0, double ry0, double rz0, double bad0) :
    badness(bad0)
{
    r[0] = rx0;
    r[1] = ry0;
    r[2] = rz0;
    badness_sum = badness;
    age = 1;
}

Atom_t::Atom_t(double r0[3], double bad0) :
    badness(bad0)
{
    copy(r0, r0+3, r);
    badness_sum = badness;
    age = 1;
}

double Atom_t::Badness() const
{
    return badness;
}

double Atom_t::AvgBadness() const
{
    return (age != 0) ? 1.0*badness_sum/age : 0.0;
}

double Atom_t::IncBadness(double db)
{
    badness += db;
    if (badness < BGA::eps_badness)
	badness = 0.0;
    badness_sum += badness;
    age++;
    return badness;
}

double Atom_t::DecBadness(double db)
{
    badness -= db;
    if (badness < BGA::eps_badness)
	badness = 0.0;
    badness_sum += badness;
    age++;
    return badness;
}

double Atom_t::ResetBadness(double b)
{
    badness = badness_sum = b;
    age = 1;
    return badness;
}

bool operator==(const Atom_t& a1, const Atom_t& a2)
{
    return equal(a1.r, a1.r+3, a2.r);
}

double dist2(const Atom_t& a1, const Atom_t& a2)
{
    BGA::cnt.distance_calls++;
    double dr[3] = { a1.r[0]-a2.r[0], a1.r[1]-a2.r[1], a1.r[2]-a2.r[2] };
    return dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
}

////////////////////////////////////////////////////////////////////////
// AtomFilter_t - definitions of subclasses
////////////////////////////////////////////////////////////////////////

bool BondAngleFilter_t::Check(Atom_t* pa, Molecule* pm, double* adist)
{
    return false;
}

////////////////////////////////////////////////////////////////////////
// PairDistance_t definitions
////////////////////////////////////////////////////////////////////////

void PairDistance_t::LockTo(Molecule* pM, Atom_t* pa1, Atom_t* pa2)
{
    double d = dist(*pa1, *pa2);
    vector<double>::iterator dnear = pM->dTarget.find_nearest(d);
    double dd = *dnear - d;
    double badness = pM->penalty(dd);
    if (badness < BGA::eps_badness)
	badness = 0.0;
    if (fabs(dd) < pM->tol_dd)
    {
	dUsed = +1.0 * (*dnear);
	pM->dTarget.erase(dnear);
    }
    else
    {
	dUsed = -1.0 * (*dnear);
    }
    double badnesshalf = badness/2.0;
    pa1->IncBadness(badnesshalf);
    pa2->IncBadness(badnesshalf);
    pM->badness += badness;
}

void PairDistance_t::Release(Molecule* pM, Atom_t* pa1, Atom_t* pa2)
{
    double badness = Badness(pM, pa1, pa2);
    double badnesshalf = badness/2.0;
    pa1->DecBadness(badnesshalf);
    pa2->DecBadness(badnesshalf);
    pM->badness -= badness;
    if (pM->badness < BGA::eps_badness)
	pM->badness = 0.0;
    if (dUsed > 0.0)
	pM->dTarget.return_back(dUsed);
}

double PairDistance_t::Badness(Molecule *pM, Atom_t *pa1, Atom_t *pa2)
{
    return pM->penalty( fabs(dUsed) - dist(*pa1, *pa2) );
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
    *this = v;
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

vector<double>::iterator DistanceTable::find_nearest(const double& dfind)
{
    iterator ii = lower_bound(begin(), end(), dfind);
    if (    ( ii == end() && size() != 0 ) ||
	    ( ii != begin() && (dfind - *(ii-1)) < (*ii - dfind) )
       )
	--ii;
    return ii;
}

vector<double>::iterator DistanceTable::return_back(const double& dback)
{
    iterator ii = lower_bound(begin(), end(), dback);
    return insert(ii, dback);
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
bool   Molecule::evolve_jump = true;
bool   Molecule::evolve_relax = false;
bool   Molecule::degenerate_relax = false;
int    Molecule::center_size = 40;
vector<AtomFilter_t> Molecule::atom_filters;

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
    for (int i = 0; i < vx.size(); ++i)
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
    atoms.resize(M.atoms.size());
    vector<Atom_t*>::const_iterator asrc;
    vector<Atom_t*>::iterator adup;
    for (asrc = M.atoms.begin(), adup = atoms.begin();
	    asrc != M.atoms.end(); ++asrc, ++adup)
    {
	*adup = new Atom_t(**asrc);
    }
    // finished duplication
    val_max_NAtoms = M.val_max_NAtoms;
    // map source atom pointers to this atom pointers
    map<const Atom_t*, Atom_t*> patom_clone;
    for (asrc = M.atoms.begin(), adup = atoms.begin();
	    asrc != M.atoms.end();  ++asrc, ++adup)
    {
	patom_clone[*asrc] = *adup;
    }
    typedef map<OrderedPair<Atom_t*>,PairDistance_t>::const_iterator MAPcit;
    for (MAPcit ii = M.pairs.begin(); ii != M.pairs.end(); ++ii)
    {
	Atom_t* a1 = patom_clone[ii->first.first];
	Atom_t* a2 = patom_clone[ii->first.second];
	OrderedPair<Atom_t*> key(a1, a2);
	pairs[key] = ii->second;
    }
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
    val_max_NAtoms = dTarget.NAtoms;
    // default output format
    OutFmtXYZ();
}

Molecule::~Molecule()
{
    // we must call Clear() to delete Atom_t objects
    Clear();
}


//////////////////////////////////////////////////////////////////////////
//// Molecule badness/fitness evaluation
//////////////////////////////////////////////////////////////////////////

double Molecule::penalty(double dd)
{
    BGA::cnt.penalty_calls++;
    return dd*dd;
}

namespace MoleculeRecalculate
{
    typedef map<OrderedPair<Atom_t*>,PairDistance_t>::iterator MAPit;
    struct BadnessWithMAPit_t
    {
	double badness;
	MAPit  iter;
    };
    bool comp_BadnessWithMAPit_t_Badness(const BadnessWithMAPit_t& lhs,
	    const BadnessWithMAPit_t& rhs)
    {
	return lhs.badness < rhs.badness;
    }
}

void Molecule::Recalculate()
{
    using namespace MoleculeRecalculate;
    if (NAtoms() > max_NAtoms())
    {
	cerr << "E: molecule too large in Recalculate()" << endl;
	throw InvalidMolecule();
    }
    // reset molecule
    badness = 0;
    // reset all atoms
    typedef vector<Atom_t*>::iterator VPAit;
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai)
    {
	(*pai)->ResetBadness();
    }
    // order pair iterators by corresponding badness for accurate summation
    // first create and fill array of BadnessWithMAPit_t
    BadnessWithMAPit_t ordered_bwi[pairs.size()];
    BadnessWithMAPit_t* pbwi = ordered_bwi;
    for (MAPit ii = pairs.begin(); ii != pairs.end(); ++ii, ++pbwi)
    {
	Atom_t* pa1 = ii->first.first;
	Atom_t* pa2 = ii->first.second;
	PairDistance_t& pd = ii->second;
	pbwi->badness = pd.Badness(this, pa1, pa2);
	pbwi->iter = ii;
    }
    // sort by badness
    sort(ordered_bwi, ordered_bwi + pairs.size(),
	    comp_BadnessWithMAPit_t_Badness);
    // sum over sorted iterators
    for (pbwi = ordered_bwi; pbwi != ordered_bwi + pairs.size(); ++pbwi)
    {
	Atom_t* pa1 = pbwi->iter->first.first;
	Atom_t* pa2 = pbwi->iter->first.second;
	double badnesshalf = pbwi->badness/2.0;
	badness += pbwi->badness;
	pa1->IncBadness(badnesshalf);
	pa2->IncBadness(badnesshalf);
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


//////////////////////////////////////////////////////////////////////////
// Molecule operators
//////////////////////////////////////////////////////////////////////////

bool operator==(const Molecule& m1, const Molecule& m2)
{
    if (&m1 == &m2)
	return true;
    bool isequal =
	m1.max_NAtoms() == m2.max_NAtoms() &&
	m1.atoms == m2.atoms;
    return isequal;
}

void Molecule::Set_max_NAtoms(int s)
{
    if (s > dTarget.NAtoms && tol_dd > 0.0)
    {
	cerr << "E: not enough distances for max_NAtoms = " << s << '.' <<
	    "  Did you forget tol_dd = 0?" << endl;
	throw InvalidMolecule();
    }
    else if (s < 1)
    {
	cerr << "E: invalid value of max_NAtoms = " << s << endl;
	throw InvalidMolecule();
    }
    else if (s < NAtoms())
    {
	cerr << "E: molecule too large in Set_max_NAtoms()" << endl;
	throw InvalidMolecule();
    }
    val_max_NAtoms = s;
}

Molecule& Molecule::Shift(double dx, double dy, double dz)
{
    typedef vector<Atom_t*>::iterator VPAit;
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai)
    {
	(*pai)->r[0] += dx;
	(*pai)->r[1] += dy;
	(*pai)->r[2] += dz;
    }
    return *this;
}

Molecule& Molecule::Center()
{
    double avg_rx = 0.0, avg_ry = 0.0, avg_rz = 0.0;
    typedef vector<Atom_t*>::iterator VPAit;
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai)
    {
	avg_rx += (*pai)->r[0];
	avg_ry += (*pai)->r[1];
	avg_rz += (*pai)->r[2];
    }
    avg_rx /= NAtoms();
    avg_ry /= NAtoms();
    avg_rz /= NAtoms();
    Shift(-avg_rx, -avg_ry, -avg_rz);
    return *this;
}

Molecule& Molecule::Pop(const int cidx)
{
    if (cidx < 0 || cidx >= NAtoms())
    {
	throw range_error("in Molecule::Pop(const int cidx)");
    }
    typedef map<OrderedPair<Atom_t*>,PairDistance_t>::iterator MAPit;
    typedef vector<Atom_t*>::iterator VPAit;
    VPAit popped = atoms.begin() + cidx;
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai)
    {
	if (pai == popped)  continue;
	OrderedPair<Atom_t*> key(*pai, *popped);
	MAPit pit = pairs.find(key);
	PairDistance_t& pd = pit->second;
	pd.Release(this, *pai, *popped);
	pairs.erase(pit);
    }
    delete *popped;
    atoms.erase(popped);
    return *this;
}

Molecule& Molecule::Pop(const list<int>& cidx)
{
    list<int> atoms2pop;
    // build a list of indices of atoms to be popped
    for ( list<int>::const_iterator ii = cidx.begin();
	    ii != cidx.end(); ++ii )
    {
	if (*ii < 0 || *ii >= NAtoms())
	{
	    cerr << "index out of range in Molecule::Pop(list<int>)" << endl;
	    throw range_error("in Molecule::Pop(list<int>)");
	}
	atoms2pop.push_back(*ii);
    }
    // now pop those atoms going from the highest index down
    atoms2pop.sort();
    atoms2pop.unique();
    for ( list<int>::reverse_iterator rii = atoms2pop.rbegin();
	    rii != atoms2pop.rend(); ++rii )
    {
	Pop(*rii);
    }
    return *this;
}

Molecule& Molecule::Clear()
{
    pairs.clear();
    typedef vector<Atom_t*>::iterator VPAit;
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai)
    {
	delete *pai;
    }
    atoms.clear();
    badness = 0.0;
    return *this;
}

Molecule& Molecule::Add(Molecule& M)
{
    typedef vector<Atom_t*>::iterator VPAit;
    for (VPAit pai = M.atoms.begin(); pai != M.atoms.end(); ++pai)
    {
	Add(**pai);
    }
    return *this;
}

Molecule& Molecule::Add(double rx0, double ry0, double rz0)
{
    Add(Atom_t(rx0, ry0, rz0));
    return *this;
}

Molecule& Molecule::Add(Atom_t atom)
{
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: molecule too large in Add()" << endl;
	throw InvalidMolecule();
    }
    Atom_t* pnew_atom;
    pnew_atom = new Atom_t(atom);
    pnew_atom->ResetBadness();
    atoms.push_back(pnew_atom);
    // pnew_atom is at the end of the list
    typedef vector<Atom_t*>::iterator VPAit;
    for (VPAit pai = atoms.begin(); *pai != pnew_atom; ++pai)
    {
	Atom_t* pmol_atom = *pai;
	OrderedPair<Atom_t*> key(pmol_atom, pnew_atom);
	PairDistance_t& new_pair = pairs[key];
	new_pair.LockTo(this, pmol_atom, pnew_atom);
    }
    return *this;
}

Atom_t Molecule::Atom(const int cidx)
{
    if (cidx < 0 || cidx >= NAtoms())
    {
	throw range_error("in Molecule::Pop(list<int>)");
    }
    return *(atoms[cidx]);
}

double Molecule::calc_test_badness(Atom_t& ta, double hi_abad)
{
    // badness is calculated exactly only when <= hi_abad
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: Molecule too large in calc_test_badness()" << endl;
	throw InvalidMolecule();
    }
    double tbad = ta.Badness();
    map<int,bool> used;
    typedef vector<Atom_t*>::iterator VPAit;
    for (   VPAit pai = atoms.begin();
	    pai != atoms.end() && tbad <= hi_abad; ++pai )
    {
	double d = dist(ta, **pai);
	int idx = dTarget.find_nearest(d) - dTarget.begin();
	typedef map<int,bool>::iterator USEit;
	// adjust idx if it is already used
	USEit idx_used_it = used.find(idx);
	if (idx_used_it != used.end())
	{
	    int hi, lo, nidx = -1;
	    int dtsize = dTarget.size();
	    USEit hi_used_it = idx_used_it;
	    for (   hi = idx+1, ++hi_used_it;
		    hi < dtsize && hi_used_it != used.end() &&
		    hi == hi_used_it->first;  ++hi, ++hi_used_it )
	    { }
	    if (hi < dtsize)
		nidx = hi;
	    typedef map<int,bool>::reverse_iterator USErit;
	    USErit lo_used_rit(idx_used_it);
	    for (   lo = idx; lo >= 0 && lo_used_rit != used.rend() &&
		    lo == lo_used_rit->first; --lo, ++lo_used_rit )
	    { }
	    if (lo >= 0 && (nidx < 0 || d-dTarget[lo] < dTarget[nidx]-d))
		nidx = lo;
	    idx = nidx;
	}
	double dd = dTarget[idx] - d;
	tbad += penalty(dd);
	if (fabs(dd) < tol_dd)
	{
	    used[idx] = true;
	}
    }
    return tbad;
}

void Molecule::filter_good_atoms(vector<Atom_t>& vta,
	double evolve_range, double hi_abad)
{
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: Molecule too large in filter_good_atoms()" << endl;
	throw InvalidMolecule();
    }
    double lo_abad = hi_abad - evolve_range;
    typedef vector<Atom_t>::iterator VAit;
    typedef vector<Atom_t*>::iterator VPAit;
    // obtain badness of every test atom and adjust hi_abad cutoff
    // badness is exact only when <= hi_abad
    for (VAit tai = vta.begin(); tai != vta.end(); ++tai)
    {
	// initial atom badness is the badness sum of base atoms
	// in the next round it is 0; in any case we can use IncBadness
	double tbad = calc_test_badness(*tai, hi_abad);
	tai->IncBadness(tbad);
	if (tai->Badness() < lo_abad)
	{
	    lo_abad = tai->Badness();
	    hi_abad = lo_abad + evolve_range;
	}
    }
    // hi_abad has a correct value here
    // let us keep only good atoms
    VAit gai = vta.begin();
    for (VAit tai = vta.begin(); tai != vta.end(); ++tai)
    {
	if (tai->Badness() <= hi_abad)
	    *(gai++) = *tai;
    }
    vta.erase(gai, vta.end());
}

//pj:	isgood = abad_filter.Check(&tai);

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
    ta.r[0] = gsl_vector_get(x, 0);
    ta.r[1] = gsl_vector_get(x, 1);
    ta.r[2] = gsl_vector_get(x, 2);
    typedef vector<Atom_t*> VPA;
    typedef valarray<double> VAD;
    VPA& atoms = *( ((struct rxa_par *)params)->atoms );
    VAD& ad0 = *( ((struct rxa_par *)params)->ad0 );
    VAD& wt = *( ((struct rxa_par *)params)->wt );
    VAD ad(ad0.size());
    typedef vector<Atom_t*>::iterator VPAit;
    int idx = 0;
    gsl_matrix_set_zero(J);
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai, ++idx)
    {
	ad[idx] = dist(ta, **pai);
	double fi = wt[idx]*(ad[idx] - ad0[idx]);
	gsl_vector_set(f, idx, fi);
	for (int jdx = 0; ad[idx] != 0.0 && jdx < 3; ++jdx)
	    gsl_matrix_set(J, idx, jdx,
		    wt[idx]*(ta.r[jdx] - (*pai)->r[jdx])/ad[idx]);
    }
    return GSL_SUCCESS;
}

int rxa_f(const gsl_vector* x, void* params, gsl_vector* f)
{
    static Atom_t ta(0.0, 0.0, 0.0);
    ta.r[0] = gsl_vector_get(x, 0);
    ta.r[1] = gsl_vector_get(x, 1);
    ta.r[2] = gsl_vector_get(x, 2);
    typedef vector<Atom_t*> VPA;
    typedef valarray<double> VAD;
    VPA& atoms = *( ((struct rxa_par *)params)->atoms );
    VAD& ad0 = *( ((struct rxa_par *)params)->ad0 );
    VAD& wt = *( ((struct rxa_par *)params)->wt );
    typedef vector<Atom_t*>::iterator VPAit;
    int idx = 0;
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai, ++idx)
    {
	double fi = wt[idx]*(dist(ta, **pai) - ad0[idx]);
	gsl_vector_set(f, idx, fi);
    }
    return GSL_SUCCESS;
}

int rxa_df(const gsl_vector* x, void* params, gsl_matrix* J)
{
    typedef vector<Atom_t*> VPA;
    typedef valarray<double> VAD;
    VPA& atoms = *( ((struct rxa_par *)params)->atoms );
    // just use rxa_fdf ignoring function values
    gsl_vector* fignore = gsl_vector_alloc(atoms.size());
    int status = rxa_fdf(x, params, fignore, J);
    gsl_vector_free(fignore);
    return status;
}

Molecule& Molecule::RelaxAtom(vector<Atom_t*>::iterator pai)
{
    RelaxAtom(pai - atoms.begin());
    return *this;
}

Molecule& Molecule::RelaxAtom(const int cidx)
{
    if (cidx < 0 || cidx >= NAtoms())
    {
	throw range_error("in Molecule::RelaxAtom(list<int>)");
    }
    Atom_t ta = *(atoms[cidx]);
    Pop(cidx);
    RelaxExternalAtom(ta);
    Add(ta);
    return *this;
}

void Molecule::RelaxExternalAtom(Atom_t& ta)
{
    const int maximum_relaxations = 20;
    // pj: this seems to be crashing when NAtoms() < 3
    if (NAtoms() < 3)
	return;
    const int max_iter = 500;
    // find if it is a member atom:
    // dTarget is not changed in this function
    int dTsize = dTarget.size();
    bool dUsed[dTsize];
    fill(dUsed, dUsed+dTsize, false);
    // prepare valarrays for rxa_* functions
    valarray<double> wt(1.0, NAtoms());
    double* pd = &wt[0];
    typedef vector<Atom_t*>::iterator VPAit;
    // first fill the array with badness
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai, ++pd)
	*pd = (*pai)->Badness();
    // convert to square root of fitness
    wt = sqrt(vdrecipw0(wt));
    // array for target distances
    valarray<double> ad0(atoms.size());
    // do relaxation on a copy of ta
    Atom_t rta(ta);
    // loop while badness is improved
    double lo_abad = numeric_limits<double>().max();
    for (int nrelax = 0; nrelax < maximum_relaxations; ++nrelax)
    {
	list<int> dUsedIdx;
	double* pad0 = &ad0[0];
	double tbad = 0.0;
	typedef vector<Atom_t*>::iterator VPAit;
	for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai)
	{
	    double d = dist(**pai, rta);
	    int idx = dTarget.find_nearest(d) - dTarget.begin();
	    if (dUsed[idx])
	    {
		int hi, lo, nidx = -1;
		for (hi = idx; hi != dTsize && dUsed[hi]; ++hi) { };
		if (hi != dTsize)
		    nidx = hi;
		for (lo = idx; lo >= 0 && dUsed[lo]; --lo) { };
		if (lo >= 0  && (nidx < 0 || d-dTarget[lo] < dTarget[nidx]-d))
		    nidx = lo;
		idx = nidx;
	    }
	    *(pad0++) = dTarget[idx];
	    double dd = dTarget[idx] - d;
	    tbad += penalty(dd);
	    if (fabs(dd) < tol_dd)
	    {
		dUsed[idx] = true;
		dUsedIdx.push_back(idx);
	    }
	}
	// restore dUsed
	for (   list<int>::iterator ii = dUsedIdx.begin();
		ii != dUsedIdx.end(); ++ii  )
	{
	    dUsed[*ii] = false;
	}
	// get out if lo_abad did not improve
	if ( eps_lt(tbad, lo_abad) )
	{
	    lo_abad = tbad;
	    ta = rta;
	    if (lo_abad < BGA::eps_badness)
		break;
	}
	else
	    break;
//	cout << "original tbad = " << tbad << endl;
	// parameter pair for minimizer functions
	struct rxa_par rxap = { &atoms, &ad0, &wt };
	// define function to be minimized
	gsl_multifit_function_fdf f;
	f.f = &rxa_f;
	f.df = &rxa_df;
	f.fdf = &rxa_fdf;
	f.n = NAtoms();
	f.p = 3;
	f.params = &rxap;
	// bind rta coordinates to vector x
	gsl_vector_view x = gsl_vector_view_array(rta.r, 3);
	// allocate solver
	gsl_multifit_fdfsolver* lms = gsl_multifit_fdfsolver_alloc(
		gsl_multifit_fdfsolver_lmsder, NAtoms(), 3);
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
	while (status == GSL_CONTINUE && iter < max_iter);
	for (int i = 0; i < 3; ++i)
	    rta.r[i] = gsl_vector_get(lms->x, i);
	// get relaxed value of tbad
	// pj: this may go away later
	/*
	tbad = 0.0;
	for (int i = 0; i < NAtoms(); ++i)
	{
	    tbad += penalty(gsl_vector_get(lms->f, i));
	}
	// and update lo_abad before reassigning the distances
	lo_abad = min(lo_abad, tbad);
	cout << "relaxed tbad = " << tbad << endl;
	*/
	gsl_vector_free(G);
	gsl_multifit_fdfsolver_free(lms);
    }
}

int Molecule::push_good_distances(
	vector<Atom_t>& vta, double* afit, int ntrials
	)
{
    // add new atom in direction defined by 2 atoms
    if (NAtoms() == max_NAtoms())
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
    // ntrials
    // N*(N-1) possible directions, Nuniqd lengths
    // check for int overflow
    double xmax_ntrials = min(
	    1.0*NAtoms()*(NAtoms()-1)*dTarget.Nuniqd + 1.0,
	    double( numeric_limits<int>().max() )
	    );
    int max_ntrials = int(xmax_ntrials);
    ntrials = min(ntrials, max_ntrials);
    int push_count = 0;
    for (int nt = 0; nt < ntrials; ++nt)
    {
	vector<int> aidx = random_wt_choose(min(NAtoms(),2), afit, NAtoms());
	Atom_t& a0 = *atoms[aidx[0]];
	valarray<double> rdir(0.0, 3);
	if (NAtoms() > 1)
	{
	    Atom_t& a1 = *atoms[aidx[1]];
	    for (int i = 0; i < 3; ++i)
		rdir[i] = a1.r[i] - a0.r[i];
	}
	// normalize rdir if defined
	bool lattice_rdir;
	double nm_rdir = vdnorm(rdir);
	if (nm_rdir > eps_d)
	{
	    rdir /= nm_rdir;
	    lattice_rdir = true;
	}
	// otherwise generate random direction
	else
	{
	    // pick a random direction
	    double phi = 2*M_PI*gsl_rng_uniform(BGA::rng);
	    double z = 2*gsl_rng_uniform(BGA::rng) - 1.0;
	    double w = sqrt(1 - z*z);
	    rdir[0] = w*cos(phi);
	    rdir[1] = w*sin(phi);
	    rdir[2] = z;
	    lattice_rdir = false;
	}
	// pick free distance
	int didx = gsl_rng_uniform_int(BGA::rng, dTarget.size());
	double radius = dTarget[didx];
	// add front atom
	double nr[3];
	for (int i = 0; i < 3; ++i)
	    nr[i] = a0.r[i] + rdir[i]*radius;
	Atom_t ad1front(nr[0], nr[1], nr[2]);
	vta.push_back(ad1front);
	++push_count;
	// check opposite direction when it makes sense
	// this accounts for extra trial
	if (lattice_rdir)
	{
	    ++nt;
	    for (int i = 0; i < 3; ++i)
		nr[i] = a0.r[i] - rdir[i]*radius;
	    Atom_t ad1back(nr[0], nr[1], nr[2]);
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
    // generate randomly oriented triangles
    if (NAtoms() == max_NAtoms())
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
    // ntrials
    // (N over 3)*4*Nuniqd^2 plane orientations and possible triangles
    // check for int overflow
    double xmax_ntrials = min(
	    4.0*(NAtoms()*(NAtoms()-1)*(NAtoms()-2)/6 *
		dTarget.Nuniqd*dTarget.Nuniqd + 2),
	    double( numeric_limits<int>().max() )
	    );
    int max_ntrials = int(xmax_ntrials);
    ntrials = min(ntrials, max_ntrials);
    int push_count = 0;
    for (int nt = 0; nt < ntrials; ++nt)
    {
	// pick 2 atoms for base and 3rd for plane orientation
	int nchoose = NAtoms() > 2 ? 3 : 2;
	vector<int> aidx = random_wt_choose(nchoose, afit, NAtoms());
	Atom_t& a0 = *atoms[aidx[0]];
	Atom_t& a1 = *atoms[aidx[1]];
	// pick 2 vertex distances
	bool with_repeat = (tol_dd <= 0.0);
	vector<int> didx = random_choose_few(2, dTarget.size(), with_repeat);
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
	valarray<double> longdir(3);
	for (int i = 0; i < 3; ++i)
	    longdir[i] = (a1.r[i] - a0.r[i])/r12;
	// generate direction perpendicular to longdir
	valarray<double> perpdir(0.0, 3);
	if (nchoose > 2)
	{
	    Atom_t& a2 = *atoms[aidx[2]];
	    for (int i = 0; i < 3; ++i)
		perpdir[i] = a2.r[i] - a0.r[i];
	    perpdir -= longdir*vddot(longdir, perpdir);
	}
	// normalize perpdir if defined
	bool lattice_plane;
	double nm_perpdir = vdnorm(perpdir);
	if (nm_perpdir > eps_d)
	{
	    perpdir /= nm_perpdir;
	    lattice_plane = true;
	}
	// otherwise generate random direction
	else
	{
	    valarray<double> pdir1(3), pdir2(3);
	    pdir1[0] = -longdir[1];
	    pdir1[1] = longdir[0];
	    pdir1[2] = 0.0;
	    if (pdir1[0] == 0 && pdir1[1] == 0)
		pdir1[0] = 1.0;
	    pdir1 /= vdnorm(pdir1);
	    pdir2 = vdcross(longdir, pdir1);
	    double phi = 2*M_PI*gsl_rng_uniform(BGA::rng);
	    perpdir = cos(phi)*pdir1 + sin(phi)*pdir2;
	    lattice_plane = false;
	}
	// allocate vallarays for positions of a0 and vertex P
	valarray<double> Pa0(a0.r, 3);
	valarray<double> P(3);
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
		Atom_t ad2(P[0], P[1], P[2]);
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
    if (NAtoms() == max_NAtoms())
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
    // ntrials
    // (N over 3)*6*2*Nuniqd^3 possible pyramids
    // check for int overflow
    double xmax_ntrials = min(
	    12.0*(NAtoms()*(NAtoms()-1)*(NAtoms()-2)/6 *
		dTarget.Nuniqd*dTarget.Nuniqd*dTarget.Nuniqd + 2),
	    double( numeric_limits<int>().max() )
	    );
    int max_ntrials = int(xmax_ntrials);
    ntrials = min(ntrials, max_ntrials);
    int push_count = 0;
    for (int nt = 0; nt < ntrials;)
    {
	// pick 3 base atoms
	vector<int> aidx = random_wt_choose(3, afit, NAtoms());
	Atom_t& a0 = *atoms[aidx[0]];
	Atom_t& a1 = *atoms[aidx[1]];
	Atom_t& a2 = *atoms[aidx[2]];
	// pick 3 vertex distances
	bool with_repeat = (tol_dd <= 0.0);
	vector<int> didx = random_choose_few(3, dTarget.size(), with_repeat);
	// loop over all permutations of selected distances
	sort(didx.begin(), didx.end());
	do
	{
	    ++nt;
	    double r14 = dTarget[ didx[0] ];
	    double r24 = dTarget[ didx[1] ];
	    double r34 = dTarget[ didx[2] ];
	    // uvi is a unit vector in a0a1 direction
	    double uvi_val[3] = { a1.r[0]-a0.r[0], a1.r[1]-a0.r[1], a1.r[2]-a0.r[2] };
	    valarray<double> uvi(uvi_val, 3);
	    double r12 = vdnorm(uvi);
	    if (r12 < eps_d)
		continue;
	    uvi /= r12;
	    // v13 is a0a2 vector
	    double v13_val[3] = { a2.r[0]-a0.r[0], a2.r[1]-a0.r[1], a2.r[2]-a0.r[2] };
	    valarray<double> v13(v13_val, 3);
	    // uvj lies in a0a1a2 plane and is perpendicular to uvi
	    valarray<double> uvj = v13;
	    uvj -= uvi*vddot(uvi, uvj);
	    double nm_uvj = vdnorm(uvj);
	    if (nm_uvj < eps_d)
		continue;
	    uvj /= nm_uvj;
	    // uvk is a unit vector perpendicular to a0a1a2 plane
	    valarray<double> uvk = vdcross(uvi, uvj);
	    double xP1 = -0.5/(r12)*(r12*r12 + r14*r14 - r24*r24);
	    // Pn are coordinates in pyramid coordinate system
	    valarray<double> P1(3);
	    P1[0] = xP1; P1[1] = P1[2] = 0.0;
	    // vT is translation from pyramid to normal system
	    valarray<double> vT(3);
	    vT[0] = a0.r[0] - xP1*uvi[0];
	    vT[1] = a0.r[1] - xP1*uvi[1];
	    vT[2] = a0.r[2] - xP1*uvi[2];
	    // obtain coordinates of P3
	    valarray<double> P3 = P1;
	    P3[0] += vddot(uvi, v13);
	    P3[1] += vddot(uvj, v13);
	    P3[2] = 0.0;
	    double xP3 = P3[0];
	    double yP3 = P3[1];
	    // find pyramid vertices
	    valarray<double> P4(0.0, 3);
	    double h2 = r14*r14 - xP1*xP1;
	    // does P4 belong to a0a1 line?
	    if (fabs(h2) < eps_d)
	    {
		// is vertex on a0a1
		if (fabs(vdnorm(P3) - r14) > eps_d)
		    continue;
		P4 = vT;
		Atom_t ad3(P4[0], P4[1], P4[2]);
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
	    vta.push_back(ad3top);
	    ++push_count;
	    // and bottom one, which makes an extra trial
	    ++nt;
	    P4 = yP4*uvj - zP4*uvk + vT;
	    Atom_t ad3bottom(P4[0], P4[1], P4[2]);
	    vta.push_back(ad3bottom);
	    push_count++;
	} while ( next_permutation(didx.begin(), didx.end()) );
    }
    return push_count;
}

Molecule& Molecule::Evolve(int ntd1, int ntd2, int ntd3)
{
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: full-sized molecule cannot Evolve()" << endl;
	throw InvalidMolecule();
    }
    vector<Atom_t> vta;
    // evolution is trivial for empty or 1-atom molecule
    switch (NAtoms())
    {
	case 0:
	    Add(0, 0, 0);
	    return *this;
	case 1:
	    double afit1 = 1.0;
	    push_good_distances(vta, &afit1, 1);
	    Add(vta[0]);
	    Center();
	    return *this;
    }
    // otherwise we need to build array of atom fitnesses
    valarray<double> vafit(NAtoms());
    double* pd = &vafit[0];
    // first fill the array with badness
    typedef vector<Atom_t*>::iterator VPAit;
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai, ++pd)
	*pd = (*pai)->Badness();
    // then get the reciprocal value
    vafit = vdrecipw0(vafit);
    double* afit = &vafit[0];
    // and push appropriate numbers of test atoms
    // here NAtoms() >= 2
    push_good_distances(vta, afit, ntd1);
    push_good_triangles(vta, afit, ntd2);
    if (NAtoms() > 2)  push_good_pyramids(vta, afit, ntd3);
    if (!vta.size())   return *this;
    // set badness range from min_badness
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
	// first fill the array with badness
	double* pfit = &vtafit[0];
	for (VAit ai = vta.begin(); ai != vta.end(); ++ai, ++pfit)
	    *pfit = ai->Badness();
	// then get the reciprocal value
	vtafit = vdrecipw0(vtafit);
	int idx = random_wt_choose(1, &vtafit[0], vtafit.size()).front();
	Add(vta[idx]);
	hi_abad = vta[idx].Badness() + evolve_range;
	vta.erase(vta.begin()+idx);
	if (evolve_relax)
	{
	    VPAit worst = max_element(atoms.begin(), atoms.end(),
		    comp_pAtom_Badness);
	    if (eps_gt((*worst)->Badness(), 0.0))
	    {
		RelaxAtom(worst);
	    }
	}
	if (NAtoms() == max_NAtoms() || !evolve_jump)
	    break;
	for (VAit ai = vta.begin(); ai != vta.end(); ++ai)
	    ai->ResetBadness();
    }
    if (NAtoms() < center_size)
	Center();
    return *this;
}

Molecule& Molecule::Degenerate(int Npop)
{
    Npop = min(NAtoms(), Npop);
    if (Npop == 0)  return *this;
    // build array of atom badnesses
    double abad[NAtoms()];
    double* pb = abad;
    typedef vector<Atom_t*>::iterator VPAit;
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai, ++pb)
    {
	*pb = (*pai)->Badness();
    }
    // generate list of atoms to pop
    vector<int> ipop_vector = random_wt_choose(Npop, abad, NAtoms());
    list<int> ipop(ipop_vector.begin(), ipop_vector.end());
    Pop(ipop);
    if (degenerate_relax && NAtoms() > 1)
    {
	VPAit worst = max_element(atoms.begin(), atoms.end(),
		comp_pAtom_Badness);
	if (eps_gt((*worst)->Badness(), 0.0))
	{
	    RelaxAtom(worst);
	}
    }
    if (NAtoms() < center_size)
	Center();
    return *this;
}

vector<int> random_choose_few(int K, int Np, bool with_repeat)
{
    vector<int> vec(K);
    if (with_repeat)
    {
	for (int i = 0; i < K; ++i)
	    vec[i] = gsl_rng_uniform_int(BGA::rng, Np);
	return vec;
    }
    // no repeats allowed here
    if (K > Np)
    {
	cerr << "random_wt_choose(): too many items to choose" << endl;
	throw range_error("too many items to choose");
    }
    // check trivial case
    else if (K == 0)
    {
	return vec;
    }
    int N = Np;
    map<int,int> tr;
    vector<int>::iterator vecit = vec.begin();
    for (int i = 0; i < K; ++i, --N)
    {
	int k = gsl_rng_uniform_int(BGA::rng, N);
	for (int c = 0;  tr.count(k) == 1;  k = tr[k], ++c)
	{
	    if (!(c < Np))
	    {
		cerr << "random_choose_few(): too many translations" << endl;
		throw runtime_error("too many translations");
	    }
	}
	*(vecit++) = k;
	tr[k] = N-1;
    }
    return vec;
}

vector<int> random_wt_choose(int K, const double* p, int Np)
{
    vector<int> vec(K);
    if (K > Np)
    {
	cerr << "random_wt_choose(): too many items to choose" << endl;
	throw range_error("too many items to choose");
    }
    // check trivial case
    else if (K == 0)
    {
	return vec;
    }
    if ( p+Np != find_if(p, p+Np, bind2nd(less<double>(),0.0)) )
    {
	cerr << "random_wt_choose(): negative choice probability" << endl;
	throw runtime_error("negative choice probability");
    }
    // now we need to do some real work
    double prob[Np];
    copy(p, p+Np, prob);
    // integer encoding
    int val[Np];
    for (int i = 0; i != Np; ++i)  val[i] = i;
    // cumulative probability
    double cumprob[Np];
    // main loop
    vector<int>::iterator vecit = vec.begin();
    for (int i = 0, Nprob = Np; i != K; ++i, --Nprob)
    {
	// calculate cumulative probability
	partial_sum(prob, prob+Nprob, cumprob);
	// if all probabilities are 0.0, set them to equal value
	if (cumprob[Nprob-1] == 0.0)
	{
	    for (int j = 0; j != Nprob; ++j)
	    {
		prob[j] = 1.0;
		cumprob[j] = (j+1.0)/Nprob;
	    }
	}
	// otherwise we can normalize cumprob
	else
	{
	    for (double* pcp = cumprob; pcp != cumprob+Nprob; ++pcp)
	    {
		*pcp /= cumprob[Nprob-1];
	    }
	}
	// now let's do binary search on cumprob:
	double r = gsl_rng_uniform(BGA::rng);
	double* pcp = upper_bound(cumprob, cumprob+Nprob, r);
	int idx = pcp - cumprob;
	*(vecit++) = val[idx];
	// overwrite this element with the last number
	prob[idx] = prob[Nprob-1];
	val[idx] = val[Nprob-1];
    }
    return vec;
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
    for (int i = 0; i < vxyz.size(); i += 3)
    {
	Add(Atom_t(vxyz[i+0], vxyz[i+1], vxyz[i+2]));
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

Molecule& Molecule::OutFmtXYZ()
{
    output_format = XYZ;
    return *this;
}

Molecule& Molecule::OutFmtAtomEye()
{
    output_format = ATOMEYE;
    return *this;
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
	    fid << "# BGA molecule format = xyz" << endl;
	    fid << "# NAtoms = " << M.NAtoms() << endl;
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
		    for (double* pr = (*pai)->r; pr != (*pai)->r + 2; ++pr)
		    {
			xyz_extremes.push_back(*pr * scale);
		    }
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

void Molecule::PrintBadness()
{
    cout << "ABadness() =";
    double mab = (*max_element(atoms.begin(), atoms.end(),
		    comp_pAtom_Badness)) -> Badness();
    bool marked = false;
    typedef vector<Atom_t*>::iterator VPAit;
    for (VPAit pai = atoms.begin(); pai != atoms.end(); ++pai)
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
    vafit = vdrecipw0(vafit);
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
