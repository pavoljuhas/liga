/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <stdexcept>
#include <limits>
#include <utility>
#include <map>
#include <unistd.h>
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

inline bool read_header(istream& fid)
{
    string dummy;
    return read_header(fid, dummy);
}

// read as many numbers as possible
template<class T> bool read_data(istream& fid, vector<T>& v)
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
    badness_sum += badness;
    age++;
    return badness;
}

double Atom_t::DecBadness(double db)
{
    badness -= db;
    if (fabs(badness) < BGA::eps_badness)
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
// Pair_t definitions
////////////////////////////////////////////////////////////////////////

Pair_t::Pair_t(Molecule *pM, Atom_t *a1, Atom_t *a2) :
    owner(pM), atom1(a1), atom2(a2)
{
    d = dist(*atom1, *atom2);
    vector<double>::iterator dnear = owner->dTarget.find_nearest(d);
    double dd = *dnear - d;
    badness = owner->penalty(dd);
    BGA::cnt.penalty_calls++;
    if (badness < BGA::eps_badness)
	badness = 0.0;
    if (fabs(dd) < owner->tol_dd)
    {
	dUsed = *dnear;
	owner->dTarget.erase(dnear);
    }
    else
	dUsed = -1.0;
    double badnesshalf = badness/2.0;
    atom1->IncBadness(badnesshalf);
    atom2->IncBadness(badnesshalf);
    owner->badness += badness;
}

// this constructor should be used only in Molecule operator=()
Pair_t::Pair_t(Molecule *pM, Atom_t *a1, Atom_t *a2, const Pair_t& p0) :
    owner(pM), atom1(a1), atom2(a2)
{
    if (owner == p0.owner)
	throw runtime_error("cloning Pair_t() of the same molecule");
    badness = p0.badness;
    dUsed = p0.dUsed;
}

Pair_t::Pair_t(const Pair_t& pair0)
{
    cerr << "call to Pair_t() copy constructor" << endl;
    throw runtime_error("call to Pair_t() copy constructor");
}

Pair_t& Pair_t::operator=(const Pair_t&)
{
    cerr << "call to Pair_t operator=()" << endl;
    throw runtime_error("call to Pair_t operator=()");
}

Pair_t::~Pair_t()
{
    atom1->DecBadness(badness);
    atom2->DecBadness(badness);
    owner->badness -= 2*badness;
    if (fabs(owner->badness) < BGA::eps_badness)
	owner->badness = 0.0;
    if (dUsed > 0.0)
	owner->dTarget.return_back(dUsed);
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
}

DistanceTable& DistanceTable::operator= (const DistanceTable& d0)
{
    if (this == &d0)
	return *this;
    resize(d0.size());
    copy(d0.begin(), d0.end(), begin());
    NAtoms = d0.NAtoms;
    max_d = d0.max_d;
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
    max_d = *rbegin();
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
int    Molecule::center_size = 40;
namespace BGA {
    double pow2(double x) {return x*x;}
    double well(double x) {return fabs(x)<Molecule::tol_dd ? 0.0 : 1.0;}
}
//double (*Molecule::penalty)(double) = BGA::well;
//double (*Molecule::penalty)(double) = fabs;
double (*Molecule::penalty)(double) = BGA::pow2;

Molecule::Molecule(const DistanceTable& dtab) : dTarget(dtab)
{
    init();
}

Molecule::Molecule(const DistanceTable& dtab,
	const int s, const double *px, const double *py, const double *pz
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
    atoms = M.atoms;
    val_max_NAtoms = M.val_max_NAtoms;
    // map of src atoms to this.atoms
    map<const Atom_t*, Atom_t*> pclone;
    list<Atom_t>::const_iterator src;
    list<Atom_t>::iterator clone;
    for (src = M.atoms.begin(), clone = atoms.begin();
	    src != M.atoms.end();  ++src, ++clone)
    {
	pclone[&(*src)] = &(*clone);
    }
    typedef map<OrderedPair<Atom_t*>,Pair_t*>::const_iterator MAPcit;
    for (MAPcit ii = M.pairs.begin(); ii != M.pairs.end(); ++ii)
    {
	Atom_t *a1 = pclone[ii->first.first];
	Atom_t *a2 = pclone[ii->first.second];
	OrderedPair<Atom_t*> key(a1, a2);
	pairs[key] = new Pair_t(this, a1, a2, *(ii->second));
    }
    badness = M.badness;
    // IO helpers
    output_format = M.output_format;
    opened_file = M.opened_file;
    return *this;
}

void Molecule::init()
{
    badness = 0;
    val_max_NAtoms = dTarget.NAtoms;
    // default output format
    OutFmtXYZ();
}

void Molecule::Set_max_NAtoms(int s)
{
    if (s > dTarget.NAtoms)
    {
	cerr << "E: max_NAtoms = " << s <<
	    " is too large for distance table" << endl;
	throw InvalidMolecule();
    }
    else if (s < 1)
    {
	cerr << "E: invalid value of max_NAtoms = " << s << endl;
	throw InvalidMolecule();
    }
    val_max_NAtoms = s;
}

Molecule::~Molecule()
{
    // we must call Clear() to delete Pair_t objects
    Clear();
}


//////////////////////////////////////////////////////////////////////////
//// Molecule badness/fitness evaluation
//////////////////////////////////////////////////////////////////////////

void Molecule::Recalculate()
{
    if (NAtoms() > max_NAtoms())
    {
	cerr << "E: molecule too large in Recalculate()" << endl;
	throw InvalidMolecule();
    }
    // destroy all pairs
    typedef map<OrderedPair<Atom_t*>,Pair_t*>::iterator MAPit;
    for (MAPit ii = pairs.begin(); ii != pairs.end(); ++ii)
    {
	delete ii->second;
    }
    pairs.clear();
    // molecule parameters
    badness = 0;
    // reset all atoms
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	ai->ResetBadness();
    }
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	LAit aj = ai;
	for (++aj; aj != atoms.end(); ++aj)
	{
	    OrderedPair<Atom_t*> key(&(*ai), &(*aj));
	    pairs[key] = new Pair_t(this, &(*ai), &(*aj));
	}
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

bool comp_Atom_r0(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.r[0] < rhs.r[0];
}

bool comp_Atom_r1(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.r[1] < rhs.r[1];
}

bool comp_Atom_r2(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.r[2] < rhs.r[2];
}


//////////////////////////////////////////////////////////////////////////
// Molecule operators
//////////////////////////////////////////////////////////////////////////

Molecule& Molecule::Shift(double dx, double dy, double dz)
{
    for (list<Atom_t>::iterator ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	ai->r[0] += dx;
	ai->r[1] += dy;
	ai->r[2] += dz;
    }
    return *this;
}

Molecule& Molecule::Center()
{
    double avg_rx = 0.0, avg_ry = 0.0, avg_rz = 0.0;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	avg_rx += ai->r[0];
	avg_ry += ai->r[1];
	avg_rz += ai->r[2];
    }
    avg_rx /= NAtoms();
    avg_ry /= NAtoms();
    avg_rz /= NAtoms();
    Shift(-avg_rx, -avg_ry, -avg_rz);
    return *this;
}

Molecule& Molecule::Pop(list<Atom_t>::iterator api)
{
    // delete all pairs containing *api
    typedef list<Atom_t>::iterator LAit;
    for (LAit ii = atoms.begin(); ii != atoms.end(); ++ii)
    {
	if (ii == api)  continue;
	OrderedPair<Atom_t*> key(&(*ii), &(*api));
	Pair_t *pp = pairs[key];
	delete pp;
	pairs.erase(key);
    }
    atoms.erase(api);
    return *this;
}

Molecule& Molecule::Pop(const int cidx)
{
    if (cidx < 0 || cidx >= NAtoms())
    {
	throw range_error("in Molecule::Pop(list<int>)");
    }
    list<Atom_t>::iterator apop = list_at(atoms, cidx);
    Pop(apop);
    return *this;
}

Molecule& Molecule::Pop(const list<int>& cidx)
{
    typedef list<Atom_t>::iterator LAit;
    list<LAit> atoms2pop;
    // build a list of iterators to atoms to be popped
    for ( list<int>::const_iterator ii = cidx.begin();
	    ii != cidx.end(); ++ii )
    {
	if (*ii < 0 || *ii >= NAtoms())
	{
	    cerr << "index out of range in Molecule::Pop(list<int>)" << endl;
	    throw range_error("in Molecule::Pop(list<int>)");
	}
	atoms2pop.push_back( list_at(atoms, *ii) );
    }
    // now pop those atoms
    for ( list<LAit>::iterator ii = atoms2pop.begin();
	    ii != atoms2pop.end(); ++ii )
    {
	Pop(*ii);
    }
    return *this;
}

Molecule& Molecule::Clear()
{
    // pairs must be destroyed before atoms;
    typedef map<OrderedPair<Atom_t*>,Pair_t*>::iterator MAPit;
    for (MAPit ii = pairs.begin(); ii != pairs.end(); ++ii)
    {
	delete ii->second;
    }
    pairs.clear();
    atoms.clear();
    badness = 0;
    return *this;
}

Molecule& Molecule::Add(Molecule& M)
{
    for ( list<Atom_t>::iterator ai = M.atoms.begin();
	    ai != M.atoms.end(); ++ai )
    {
	Add(*ai);
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
    list<Atom_t>::iterator this_atom, ai;
    this_atom = atoms.insert(atoms.end(), atom);
    this_atom->ResetBadness();
    // this_atom is at the end of the list
    for (ai = atoms.begin(); ai != this_atom; ++ai)
    {
	OrderedPair<Atom_t*> key(&(*ai), &(*this_atom));
	pairs[key] = new Pair_t(this, &(*ai), &(*this_atom));
    }
    return *this;
}

void Molecule::calc_test_badness(Atom_t& ta)
{
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: Molecule too large in calc_test_badness()" << endl;
	throw InvalidMolecule();
    }
    double tbad = 0;
    typedef vector<double>::iterator VDit;
    DistanceTable local_dTarget(dTarget);
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	double d = dist(*ai, ta);
	VDit dnear = local_dTarget.find_nearest(d);
	double dd = *dnear - d;
	tbad += penalty(dd);
	BGA::cnt.penalty_calls++;
	if (fabs(dd) < tol_dd)
	{
	    local_dTarget.erase(dnear);
	}
    }
    ta.ResetBadness(tbad);
}

void Molecule::filter_good_atoms(vector<Atom_t>& vta,
	double evolve_range, double lo_abad)
{
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: Molecule too large in filter_good_atoms()" << endl;
	throw InvalidMolecule();
    }
    // local copy of dTarget
    DistanceTable ldTarget(dTarget);
    int ldTsize = ldTarget.size();
    bool ldUsed[ldTsize];
    fill(ldUsed, ldUsed+ldTsize, false);
    double hi_abad = lo_abad + evolve_range;
    typedef vector<Atom_t>::iterator VAit;
    typedef list<Atom_t>::iterator LAit;
    typedef vector<double>::iterator VDit;
    // loop over all test atoms
    for (VAit ta = vta.begin(); ta != vta.end(); ++ta)
    {
	double tbad = ta->Badness();
	list<int> ldUsedIdx;
	// fast, possibly incomplete badness evaluation
	for (   LAit ma = atoms.begin();
		ma != atoms.end() && tbad <= hi_abad; ++ma )
	{
	    double d = dist(*ma, *ta);
	    int idx = ldTarget.find_nearest(d) - ldTarget.begin();
	    if (ldUsed[idx])
	    {
		int hi, lo, nidx = -1;
		for (hi = idx; hi != ldTsize && ldUsed[hi]; ++hi) { };
		if (hi != ldTsize)
		    nidx = hi;
		for (lo = idx; lo >= 0 && ldUsed[lo]; --lo) { };
		if (lo >= 0  && (nidx < 0 || d-ldTarget[lo] < ldTarget[nidx]-d))
		    nidx = lo;
		idx = nidx;
	    }
	    double dd = ldTarget[idx] - d;
	    tbad += penalty(dd);
	    BGA::cnt.penalty_calls++;
	    if (fabs(dd) < tol_dd)
	    {
		ldUsed[idx] = true;
		ldUsedIdx.push_back(idx);
	    }
	}
	ta->ResetBadness(tbad);
	if (tbad < lo_abad)
	{
	    lo_abad = tbad;
	    hi_abad = lo_abad + evolve_range;
	}
	// restore ldUsed
	for (   list<int>::iterator ii = ldUsedIdx.begin();
		ii != ldUsedIdx.end(); ++ii  )
	{
	    ldUsed[*ii] = false;
	}
    }
    // now keep only good atoms
    VAit gai = vta.begin();
    for (VAit ta = vta.begin(); ta != vta.end(); ++ta)
    {
	if (ta->Badness() <= hi_abad)
	    *(gai++) = *ta;
    }
    vta.erase(gai, vta.end());
}

int rxa_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
    static Atom_t ta(0.0, 0.0, 0.0);
    ta.r[0] = gsl_vector_get(x, 0);
    ta.r[1] = gsl_vector_get(x, 1);
    ta.r[2] = gsl_vector_get(x, 2);
    typedef list<Atom_t> LA;
    typedef valarray<double> VAD;
    pair<LA*,VAD*>& rxa_par = *( (pair<LA*,VAD*> *)params );
    LA& atoms = *(rxa_par.first);
    VAD& ad0 = *(rxa_par.second);
    VAD ad(ad0.size());
    typedef list<Atom_t>::iterator LAit;
    int idx = 0;
    gsl_matrix_set_zero(J);
    for (LAit ma = atoms.begin(); ma != atoms.end(); ++ma, ++idx)
    {
	ad[idx] = dist(ta, *ma);
	gsl_vector_set(f, idx, ad[idx]-ad0[idx]);
	for (int jdx = 0; ad[idx] != 0.0 && jdx < 3; ++jdx)
	    gsl_matrix_set(J, idx, jdx, (ta.r[jdx] - ma->r[jdx])/ad[idx]);
    }
    return GSL_SUCCESS;
}

int rxa_f(const gsl_vector *x, void *params, gsl_vector *f)
{
    static Atom_t ta(0.0, 0.0, 0.0);
    ta.r[0] = gsl_vector_get(x, 0);
    ta.r[1] = gsl_vector_get(x, 1);
    ta.r[2] = gsl_vector_get(x, 2);
    typedef list<Atom_t> LA;
    typedef valarray<double> VAD;
    pair<LA*,VAD*>& rxa_par = *( (pair<LA*,VAD*> *)params );
    LA& atoms = *(rxa_par.first);
    VAD& ad0 = *(rxa_par.second);
    typedef list<Atom_t>::iterator LAit;
    int idx = 0;
    for (LAit ma = atoms.begin(); ma != atoms.end(); ++ma, ++idx)
    {
	gsl_vector_set(f, idx, dist(ta, *ma) - ad0[idx]);
    }
    return GSL_SUCCESS;
}

int rxa_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
    typedef list<Atom_t> LA;
    typedef valarray<double> VAD;
    pair<LA*,VAD*>& rxa_par = *( (pair<LA*,VAD*> *)params );
    LA& atoms = *(rxa_par.first);
    // just use rxa_fdf ignoring function values
    gsl_vector *fignore = gsl_vector_alloc(rxa_par.first->size());
    int status = rxa_fdf(x, params, fignore, J);
    gsl_vector_free(fignore);
    return status;
}

void Molecule::relax_atom(Atom_t& ta)
{
    const int max_iter = 500;
    // dTarget is not changed in this function
    int dTsize = dTarget.size();
    bool dUsed[dTsize];
    fill(dUsed, dUsed+dTsize, false);
    valarray<double> ad0(atoms.size());
    typedef list<Atom_t>::iterator LAit;
    typedef vector<double>::iterator VDit;
    // do relaxation on a copy of ta
    Atom_t rta(ta);
    // loop while badness is improved
    double lo_abad = numeric_limits<double>().max();
    while (true)
    {
	list<int> dUsedIdx;
	double *pad0 = &ad0[0];
	double tbad = 0.0;
	for (LAit ma = atoms.begin(); ma != atoms.end(); ++ma)
	{
	    double d = dist(*ma, rta);
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
	    BGA::cnt.penalty_calls++;
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
	if (tbad < lo_abad)
	{
	    lo_abad = tbad;
	    ta = rta;
	    if (lo_abad < BGA::eps_badness)
		break;
	}
	else
	    break;
	cout << "original tbad = " << tbad << endl;
	// parameter pair for minimizer functions
	typedef list<Atom_t> LA;
	typedef valarray<double> VAD;
	pair<LA*,VAD*> rxa_par(&atoms, &ad0);
	// define function to be minimized
	gsl_multifit_function_fdf f;
	f.f = &rxa_f;
	f.df = &rxa_df;
	f.fdf = &rxa_fdf;
	f.n = NAtoms();
	f.p = 3;
	f.params = &rxa_par;
	// bind rta coordinates to vector x
	gsl_vector_view x = gsl_vector_view_array(rta.r, 3);
	// allocate solver
	gsl_multifit_fdfsolver *lms = gsl_multifit_fdfsolver_alloc(
		gsl_multifit_fdfsolver_lmsder, NAtoms(), 3);
	gsl_multifit_fdfsolver_set(lms, &f, &x.vector);
	// minimize atom badness
	int iter = 0, status;
	do
	{
	    ++iter;
	    status = gsl_multifit_fdfsolver_iterate(lms);
	    if (status)
	    {
		cerr << "LM solver status = " << gsl_strerror(status) << endl;
		break;
	    }
	    status = gsl_multifit_test_delta(lms->dx, lms->x, tol_r, tol_r);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	// get relaxed value of tbad
	// pj: this may go away later
	tbad = 0.0;
	for (int i = 0; i < NAtoms(); ++i)
	{
	    tbad += penalty(gsl_vector_get(lms->f, i));
	    BGA::cnt.penalty_calls++;
	}
	// and update lo_abad before reassigning the distances
//	lo_abad = min(lo_abad, tbad);
	cout << "relaxed tbad = " << tbad << endl;
	gsl_multifit_fdfsolver_free(lms);
    }
}

int Molecule::push_good_distances(
	vector<Atom_t>& vta, double *afit, int ntrials
	)
{
    // prepare discrete generator
    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(NAtoms(), afit);
    for (int nt = 0; nt < ntrials; ++nt)
    {
	// pick a random direction
	double phi = 2*M_PI*gsl_rng_uniform(BGA::rng);
	double z = 2*gsl_rng_uniform(BGA::rng) - 1.0;
	double w = sqrt(1 - z*z);
	double rdir[3];
	rdir[0] = w*cos(phi);
	rdir[1] = w*sin(phi);
	rdir[2] = z;
	// pick one atom and free distance
	int aidx = gsl_ran_discrete(BGA::rng, table);
	int didx = gsl_rng_uniform_int(BGA::rng, dTarget.size());
	double radius = dTarget[didx];
	Atom_t& a1 = *list_at(atoms, aidx);
	double nrx = a1.r[0] + rdir[0]*radius;
	double nry = a1.r[1] + rdir[1]*radius;
	double nrz = a1.r[2] + rdir[2]*radius;
	Atom_t ad1(nrx, nry, nrz);
	ad1.IncBadness(a1.Badness());
	vta.push_back(ad1);
    }
    gsl_ran_discrete_free(table);
    return ntrials;
}

int Molecule::push_good_triangles(
	vector<Atom_t>& vta, double *afit, int ntrials
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
    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(NAtoms(), afit);
    int push_count = 0;
    for (int nt = 0; nt < ntrials; ++nt)
    {
	// pick first atom and free distance
	list<int> aidx = random_wt_choose(2, afit, NAtoms());
	list<int>::iterator aidxit = aidx.begin();
	Atom_t& a1 = *list_at(atoms, *(aidxit++));
	Atom_t& a2 = *list_at(atoms, *(aidxit++));
	int idf1 = gsl_rng_uniform_int(BGA::rng, dTarget.size());
	int idf2 = gsl_rng_uniform_int(BGA::rng, dTarget.size()-1) + 1;
	idf2 = (idf2 + idf1) % dTarget.size();
	double r13 = dTarget[idf1];
	double r23 = dTarget[idf2];
	double r12 = dist(a1, a2);
	// is triangle base reasonably large?
	if (r12 < 1.0)
	    continue;
	double xlong = (r13*r13 + r12*r12 - r23*r23) / (2.0*r12);
	double xperp2 = r13*r13 - xlong*xlong;
	double xperp;
	if (xperp2 > 0.0)
	    xperp = sqrt(xperp2);
	else if (xperp2 > -0.25)
	    xperp = 0.0;
	else
	    continue;
	// direction along triangle base:
	double longdir[3] = {
	    (a2.r[0] - a1.r[0])/r12,
	    (a2.r[1] - a1.r[1])/r12,
	    (a2.r[2] - a1.r[2])/r12 };
	// generate random direction in the plane perpendicular to longdir
	double pdir1[3] = { -longdir[1],  longdir[0],  0.0 };
	if (pdir1[0] == 0 && pdir1[1] == 0)
	    pdir1[0] = 1.0;
	double norm_pdir1 =
	    sqrt(pdir1[0]*pdir1[0] + pdir1[1]*pdir1[1] + pdir1[2]*pdir1[2]);
	pdir1[0] /= norm_pdir1; pdir1[1] /= norm_pdir1; pdir1[2] /= norm_pdir1;
	double pdir2[3] = {
	    longdir[1]*pdir1[2]-longdir[2]*pdir1[1],
	    longdir[2]*pdir1[0]-longdir[0]*pdir1[2],
	    longdir[0]*pdir1[1]-longdir[1]*pdir1[0] };
	double phi = 2*M_PI*gsl_rng_uniform(BGA::rng);
	double perpdir[3] = {
	    cos(phi)*pdir1[0]+sin(phi)*pdir2[0],
	    cos(phi)*pdir1[1]+sin(phi)*pdir2[1],
	    cos(phi)*pdir1[2]+sin(phi)*pdir2[2] };
	// prepare new atom
	double nrx = a1.r[0] + xlong*longdir[0] + xperp*perpdir[0];
	double nry = a1.r[1] + xlong*longdir[1] + xperp*perpdir[1];
	double nrz = a1.r[2] + xlong*longdir[2] + xperp*perpdir[2];
	Atom_t ad2(nrx, nry, nrz);
	ad2.IncBadness(a1.Badness() + a2.Badness());
	vta.push_back(ad2);
	push_count++;
    }
    gsl_ran_discrete_free(table);
    return push_count;
}

int Molecule::push_good_pyramids(
	vector<Atom_t>& vta, double *afit, int ntrials
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
    // (N over 3)*6 distance permutations
    int max_ntrials = NAtoms()*(NAtoms()-1)*(NAtoms()-2);
    ntrials = min(ntrials, max_ntrials);
    int push_count = 0;
    for (int nt = 0; nt < ntrials;)
    {
	// pick 3 base atoms
	list<int> aidx = random_wt_choose(3, afit, NAtoms());
	list<int>::iterator aidxit = aidx.begin();
	Atom_t& a1 = *list_at(atoms, *(aidxit++));
	Atom_t& a2 = *list_at(atoms, *(aidxit++));
	Atom_t& a3 = *list_at(atoms, *(aidxit++));
	double base_badness = a1.Badness()+a2.Badness()+a3.Badness();
	// pick 3 vertex distances
	list<int> didx = random_choose_few(3, dTarget.size());
	list<int>::iterator didxit = didx.begin();
	double dvperm[3] = {
	    dTarget[*(didxit++)],
	    dTarget[*(didxit++)],
	    dTarget[*(didxit++)] };
	sort(dvperm, dvperm+3);
	// loop over unique permutations of dvperm
	do
	{
	    ++nt;
	    double r14 = dvperm[0];
	    double r24 = dvperm[1];
	    double r34 = dvperm[2];
	    // uvi is a unit vector in a1a2 direction
	    double uvi_val[3] = { a2.r[0]-a1.r[0], a2.r[1]-a1.r[1], a2.r[2]-a1.r[2] };
	    valarray<double> uvi(uvi_val, 3);
	    double r12 = vdnorm(uvi);
	    if (r12 < 1.0)
		continue;
	    uvi /= r12;
	    // v13 is a1a3 vector
	    double v13_val[3] = { a3.r[0]-a1.r[0], a3.r[1]-a1.r[1], a3.r[2]-a1.r[2] };
	    valarray<double> v13(v13_val, 3);
	    // uvj lies in a1a2a3 plane and is perpendicular to uvi
	    valarray<double> uvj = v13;
	    uvj -= uvi*vddot(uvi, uvj);
	    double nm_uvj = vdnorm(uvj);
	    if (nm_uvj < 1.0)
		continue;
	    uvj /= nm_uvj;
	    // uvk is a unit vector perpendicular to a1a2a3 plane
	    valarray<double> uvk = vdcross(uvi, uvj);
	    double xP1 = -0.5/(r12)*(r12*r12 + r14*r14 - r24*r24);
	    // Pn are coordinates in pyramid coordinate system
	    valarray<double> P1(3);
	    P1[0] = xP1; P1[1] = P1[2] = 0.0;
	    // vT is translation from pyramid to normal system
	    valarray<double> vT(3);
	    vT[0] = a1.r[0] - xP1*uvi[0];
	    vT[1] = a1.r[1] - xP1*uvi[1];
	    vT[2] = a1.r[2] - xP1*uvi[2];
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
	    // does P4 belong to a1a2 line?
	    if (fabs(h2) < 0.25)
	    {
		// is vertex on a1a2
		if (fabs(vdnorm(P3)) - r14 > 0.25)
		    continue;
		P4 = vT;
		Atom_t ad3(P4[0], P4[1], P4[2]);
		ad3.IncBadness(base_badness);
		vta.push_back(ad3);
		push_count++;
		continue;
	    }
	    else if (h2 < 0)
	    {
		continue;
	    }
	    double yP4 = 0.5/(yP3)*(h2 + xP3*xP3 + yP3*yP3 - r34*r34);
	    double z2P4 = h2 - yP4*yP4;
	    // does P4 belong to a1a2a3 plane?
	    if (fabs(z2P4) < 0.25)
	    {
		P4 = yP4*uvj + vT;
		Atom_t ad3(P4[0], P4[1], P4[2]);
		ad3.IncBadness(base_badness);
		vta.push_back(ad3);
		push_count++;
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
	    ad3top.IncBadness(base_badness);
	    vta.push_back(ad3top);
	    push_count++;
	    // and bottom one
	    P4 = yP4*uvj - zP4*uvk + vT;
	    Atom_t ad3bottom(P4[0], P4[1], P4[2]);
	    ad3bottom.IncBadness(base_badness);
	    vta.push_back(ad3bottom);
	    push_count++;
	} while (next_permutation(dvperm, dvperm+3));
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
    double *pd = &vafit[0];
    typedef list<Atom_t>::iterator LAit;
    // first fill the array with badness
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai, ++pd)
	*pd = ai->Badness();
    // then get the reciprocal value
    vafit = vdrecipw0(vafit);
    double *afit = &vafit[0];
    // and push appropriate numbers of test atoms
    // here NAtoms() >= 2
    push_good_distances(vta, afit, ntd1);
    push_good_triangles(vta, afit, ntd2);
    if (NAtoms() > 2)  push_good_pyramids(vta, afit, ntd3);
    if (!vta.size())   return *this;
    // set badness range from min_badness
    double evolve_range = NAtoms()*tol_nbad*evolve_frac;
    // try to add as many atoms as possible
    typedef vector<Atom_t>::iterator VAit;
    double lo_abad = numeric_limits<double>().max();
    while (true)
    {
	filter_good_atoms(vta, evolve_range, lo_abad);
	if (vta.size() == 0)
	    break;
	// calculate fitness of test atoms
	valarray<double> vtafit(vta.size());
	// first fill the array with badness
	double *pd = &vtafit[0];
	for (VAit ai = vta.begin(); ai != vta.end(); ++ai, ++pd)
	    *pd = ai->Badness();
	// then get the reciprocal value
	vtafit = vdrecipw0(vtafit);
	int idx = *(random_wt_choose(1, &vtafit[0], vtafit.size()).begin());
	Add(vta[idx]);
	lo_abad = vta[idx].Badness();
	vta.erase(vta.begin()+idx);
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
    double *pb = abad;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai, ++pb)
    {
	*pb = ai->Badness();
    }
    // generate list of atoms to pop
    list<int> ipop = random_wt_choose(Npop, abad, NAtoms());
    Pop(ipop);
    if (NAtoms() < center_size)
	Center();
    return *this;
}

list<int> random_choose_few(int K, int Np)
{
    list<int> lst;
    if (K > Np)
    {
	cerr << "random_wt_choose(): too many items to choose" << endl;
	throw range_error("too many items to choose");
    }
    // check trivial case
    else if (K == 0)
    {
	return lst;
    }
    int N = Np;
    map<int,int> tr;
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
	lst.push_back(k);
	tr[k] = N-1;
    }
    return lst;
}

list<int> random_wt_choose(int K, const double *p, int Np)
{
    list<int> lst;
    if (K > Np)
    {
	cerr << "random_wt_choose(): too many items to choose" << endl;
	throw range_error("too many items to choose");
    }
    // check trivial case
    else if (K == 0)
    {
	return lst;
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
    int Nprob = Np;
    // main loop
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
	    for (double *pcp = cumprob; pcp != cumprob+Nprob; ++pcp)
	    {
		*pcp /= cumprob[Nprob-1];
	    }
	}
	// now let's do binary search on cumprob:
	double r = gsl_rng_uniform(BGA::rng);
	double *pcp = upper_bound(cumprob, cumprob+Nprob, r);
	int idx = pcp - cumprob;
	lst.push_back(val[idx]);
	// overwrite this element with the last number
	prob[idx] = prob[Nprob-1];
	val[idx] = val[Nprob-1];
    }
    return lst;
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

template<class T> bool Molecule::ParseHeader::read_token(
	const char *token, T& value
	)
{
    const char *fieldsep = ":= ";
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

bool write_file(const char* file, Molecule& M)
{
    // test if file is writeable
    ofstream fid(file, ios_base::out|ios_base::ate );
    if (!fid)
    {
	ostringstream oss;
	oss << "write_file(): unable to write to '" << file << "'";
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
    bool result = (fid << M);
    fid.close();
    rename(writefile, file);
    return result;
}

bool Molecule::WriteXYZ(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtXYZ();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
    return result;
}

bool Molecule::WriteAtomEye(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtAtomEye();
    bool result = write_file(file, *this);
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
    bool result;
    switch (ph.format)
    {
	case M.XYZ:
	    result = M.ReadXYZ(fid);
	    break;
	case M.ATOMEYE:
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
    typedef list<Atom_t>::iterator LAit;
    LAit afirst = M.atoms.begin();
    LAit alast = M.atoms.end();
    switch (M.output_format)
    {
	case M.XYZ:
	    fid << "# BGA molecule format = xy" << endl;
	    fid << "# NAtoms = " << M.NAtoms() << endl;
	    for (LAit ai = afirst; ai != alast; ++ai)
	    {
		fid <<
		    ai->r[0] << '\t' <<
		    ai->r[1] << '\t' <<
		    ai->r[2] << endl;
	    }
	    break;
	case M.ATOMEYE:
	    double xyz_lo = 0.0;
	    double xyz_hi = 1.0;
	    double xyz_range = xyz_hi - xyz_lo;
	    if (M.NAtoms() > 0)
	    {
		const double scale = 1.01;
		double xyz_extremes[10] = {
		    -M.dTarget.max_d,
		    scale * min_element(afirst, alast, comp_Atom_r0)->r[0],
		    scale * min_element(afirst, alast, comp_Atom_r1)->r[1],
		    scale * min_element(afirst, alast, comp_Atom_r2)->r[2],
		    +M.dTarget.max_d,
		    scale * max_element(afirst, alast, comp_Atom_r0)->r[0],
		    scale * max_element(afirst, alast, comp_Atom_r1)->r[1],
		    scale * max_element(afirst, alast, comp_Atom_r2)->r[2],
		    // make atomeye happy
		    -1.75, 1.75,
		};
		xyz_lo = *min_element(xyz_extremes, xyz_extremes+10);
		xyz_hi = *max_element(xyz_extremes, xyz_extremes+10);
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
	    for (LAit ai = afirst; ai != alast; ++ai)
	    {
		fid <<
		    (ai->r[0] - xyz_med) / xyz_range + 0.5 << " " <<
		    (ai->r[1] - xyz_med) / xyz_range + 0.5 << " " <<
		    (ai->r[2] - xyz_med) / xyz_range + 0.5 << " " <<
		    ai->Badness() << endl;
	    }
	    break;
    }
    return fid;
}

void Molecule::PrintBadness()
{
    cout << "ABadness() =";
    double mab = max_element(atoms.begin(), atoms.end(),
		    comp_Atom_Badness) -> Badness();
    bool marked = false;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	cout << ' ';
	if (!marked && ai->Badness() == mab)
	{
	    cout << '+';
	    marked = true;
	}
	cout << ai->Badness();
    }
    cout << endl;
}

void Molecule::PrintFitness()
{
    valarray<double> vafit(NAtoms());
    double *pd = &vafit[0];
    typedef list<Atom_t>::iterator LAit;
    // first fill the array with badness
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai, ++pd)
	*pd = ai->Badness();
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
