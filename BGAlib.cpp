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
#include <gsl/gsl_randist.h>
#include "BGAlib.hpp"

// random number generator
gsl_rng * BGA::rng = gsl_rng_alloc(gsl_rng_default);

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

// read as many numbers as possible
template<class T> typename list<T>::iterator list_at(list<T>& lst, int n)
{
    typename list<T>::iterator ii;
    if (n <= lst.size()/2)
    {
	ii = lst.begin();
	advance(ii, n);
    }
    else
    {
	ii = lst.end();
	advance(ii, n-lst.size());
    }
    return ii;
}



////////////////////////////////////////////////////////////////////////
// SandSphere definitions
////////////////////////////////////////////////////////////////////////

SandSphere::SandSphere(int GridMax, const vector<double>& vt) :
    gridmax(GridMax)
{
    init(vt);
}

SandSphere::SandSphere(int GridMax, int s, const double *pt) :
    gridmax(GridMax)
{
    vector<double> vt(s);
    for (int i = 0; i < s; ++i) { vt[i] = pt[i]; }
    init(vt);
}

SandSphere::SandSphere(int GridMax, const char *file) :
    gridmax(GridMax)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    // read values to vt
    vector<double> vt;
    bool result = read_header(fid) && read_data(fid, vt);
    // check if everything was read
    if ( !result || !fid.eof() )
    {
	fid.clear();
	cerr << "E: " << file << ':' << fid.tellg() <<
	    ": error reading SandSphere" << endl;
	throw IOError();
    }
    fid.close();
    init(vt);
}

void SandSphere::init(const vector<double>& t)
{
    NDist = t.size();
    if (NDist == 0)
    {
	cerr << "E: target distance table is empty" << endl;
	throw InvalidDistanceTable();
    }
    // calculate and check NAtoms
    double xNAtoms = 0.5 + sqrt(1 + 8.0*NDist)/2.0;
    NAtoms = int(xNAtoms);
    if (double(NAtoms) != xNAtoms)
    {
	cerr << "E: incorrect length of target distance table, NAtoms=" <<
		xNAtoms << '\n';
	throw InvalidDistanceTable();
    }
    // fill in and check distance valarray d
    d.resize(NDist);
    copy(t.begin(), t.end(), &d[0]);
    sort(&d[0], &d[d.size()]);
    if (d[0] <= 0)
    {
	cerr << "E: non-positive entry in target distance table, " <<
	    "d[0]=" << d[0] << '\n';
	throw InvalidDistanceTable();
    }
    // calculate grid parameters
    dmax = d[d.size()-1];
    delta = dmax/gridmax;
    // calculate d2, d2lo, d2hi
    d2.resize(NDist);
    d2lo.resize(NDist);
    d2hi.resize(NDist);
    d /= delta;
    d2 = d*d;
    // take care of grid tolerance:
    SetGridTol(defGridTol);
}

void SandSphere::SetGridTol(double t)
{
    vGridTol = t;
    for (int i = 0; i < NDist; ++i)
    {
	d2lo[i] = (d[i] < vGridTol) ? 0.0 : pow(d[i] - vGridTol, 2);
	d2hi[i] = pow(d[i] + vGridTol, 2);
    }
    for ( list<Molecule *>::const_iterator i = molecules.begin();
	    i != molecules.end(); ++i )
    {
	(*i)->Recalculate();
    }
}

double SandSphere::GridTol()
{
    return vGridTol;
}


////////////////////////////////////////////////////////////////////////
// Atom_t definitions
////////////////////////////////////////////////////////////////////////

Atom_t::Atom_t(int h0, int k0, int l0, int bad0) :
    h(h0), k(k0), l(l0), badness(bad0)
{
    badness_sum = badness;
    age = 1;
}

int Atom_t::Badness() const
{
    return badness;
}

double Atom_t::AvgBadness() const
{
    return (age != 0) ? 1.0*badness_sum/age : 0.0;
}

int Atom_t::IncBadness(int b)
{
    badness += b;
    badness_sum += badness;
    age++;
    return badness;
}

int Atom_t::DecBadness(int b)
{
    badness -= b;
    badness_sum += badness;
    age++;
    return badness;
}

int Atom_t::ResetBadness()
{
    badness = badness_sum = 0;
    age = 1;
    return badness;
}

int dist2(const Atom_t& a1, const Atom_t& a2)
{
    int dr[3] = { (a1.h - a2.h), (a1.k - a2.k), (a1.l - a2.l) };
    return dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
}

////////////////////////////////////////////////////////////////////////
// Pair_t definitions
////////////////////////////////////////////////////////////////////////

bool comp_find_in_array(const int& lhs, const pair<double,double*>& rhs)
{
    const double value = rhs.first;
    const double *a = rhs.second;
    return a[lhs] < value;
}

Pair_t::Pair_t(Molecule *pM, Atom_t& a1, Atom_t& a2) :
    owner(pM), atom1(&a1), atom2(&a2)
{
    d2 = dist2(*atom1, *atom2);
    d = sqrt(1.0*d2);
    badness = 0;
    if (d < BGA::min_distance)
    {
	badness += BGA::min_distance_penalty;
    }
    // abbreviation
    list<int>::iterator ii;
    double *d2lo = &(owner->ss->d2lo[0]);
    double *d2hi = &(owner->ss->d2hi[0]);
    pair<double,double*> d2_d2hi(d2, d2hi);
    ii = lower_bound( owner->ssdIdxFree.begin(), owner->ssdIdxFree.end(),
	    d2_d2hi, comp_find_in_array );
    if (ii == owner->ssdIdxFree.end() || d2lo[*ii] > d2)
    {
	ssdIdxUsed = -1;
	badness++;
    }
    else
    {
	ssdIdxUsed = *ii;
	owner->ssdIdxFree.erase(ii);
    }
    atom1->IncBadness(badness);
    atom2->IncBadness(badness);
    owner->badness += 2*badness;
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
    int mab = max(atom1->Badness(), atom2->Badness());
    atom1->DecBadness(badness);
    atom2->DecBadness(badness);
    owner->badness -= 2*badness;
    if (ssdIdxUsed >= 0)
    {
	// return it back to owner->ssdIdxFree
	list<int>::iterator ii;
	ii = lower_bound( owner->ssdIdxFree.begin(), owner->ssdIdxFree.end(),
		ssdIdxUsed );
	owner->ssdIdxFree.insert(ii, ssdIdxUsed);
    }
}


////////////////////////////////////////////////////////////////////////
// Molecule definitions
////////////////////////////////////////////////////////////////////////

Molecule::Molecule(SandSphere *SS) : ss(SS)
{
    init();
}

Molecule::Molecule(SandSphere *SS,
	const int s, const int *ph, const int *pk, const int *pl) : ss(SS)
{
    init();
    for (int i = 0; i < s; ++i)
    {
	Add(ph[i], pk[i], pl[i]);
    }
}

Molecule::Molecule(SandSphere *SS,
	const vector<int>& vh, const vector<int>& vk, const vector<int>& vl
	) : ss(SS)
{
    init();
    if (vh.size() != vk.size() || vh.size() != vl.size())
    {
	cerr << "E: invalid coordinate vectors" << endl;
	throw InvalidMolecule();
    }
    for (int i = 0; i < vh.size(); ++i)
    {
	Add(vh[i], vk[i], vl[i]);
    }
}

Molecule::Molecule(SandSphere *SS,
	const int s, const double *px, const double *py, const double *pz
	) : ss(SS)
{
    init();
    int h, k, l;
    for (int i = 0; i < s; ++i)
    {
	h = (int) round(px[i] / ss->delta);
	k = (int) round(py[i] / ss->delta);
	l = (int) round(pz[i] / ss->delta);
	Add(h, k, l);
    }
}

Molecule::Molecule(SandSphere *SS,
	const vector<double>& vx, const vector<double>& vy,
	const vector<double>& vz
	) : ss(SS)
{
    init();
    if (vx.size() != vy.size() || vx.size() != vz.size())
    {
	cerr << "E: invalid coordinate vectors" << endl;
	throw InvalidMolecule();
    }
    int h, k, l;
    for (int i = 0; i < vx.size(); ++i)
    {
	h = (int) round(vx[i] / ss->delta);
	k = (int) round(vy[i] / ss->delta);
	l = (int) round(vz[i] / ss->delta);
	Add(h, k, l);
    }
}

Molecule::Molecule(const Molecule& M) :
    ss(M.ss)
{
    init();
    *this  = M;
}

Molecule& Molecule::operator=(const Molecule& M)
{
    if (this == &M) return *this;
    // Clear() must be the first statement
    Clear();
    if (ss != M.ss)
    {
	ss->molecules.remove(this);
	ss = M.ss;
	ss->molecules.push_back(this);
	ssdIdxFree.clear();
	for (int i = 0; i < ss->NDist; ++i)
	{
	    ssdIdxFree.push_back(i);
	}
    }
    atoms = M.atoms;
    Recalculate();
    // IO helpers
    output_format = M.output_format;
    opened_file = M.opened_file;
    return *this;
}

void Molecule::init()
{
    ss->molecules.push_back(this);
    max_NAtoms = ss->NAtoms;
    max_abad = 0;
    badness = 0;
    // prepare ssdIdxFree
    for (int i = 0; i < ss->NDist; ++i)
    {
	ssdIdxFree.push_back(i);
    }
    // default output format
    OutFmtGrid();
}

Molecule::~Molecule()
{
    // debug: cout << "ss->molecules.size() = " << ss->molecules.size() << endl;
    Clear();
    ss->molecules.remove(this);
}


//////////////////////////////////////////////////////////////////////////
//// Molecule badness/fitness evaluation
//////////////////////////////////////////////////////////////////////////
//
//double Molecule::dist(const int& i, const int& j) const
//{
////    return sqrt(1.0*dist2(i, j));
//    return 0;
//}
//
//namespace BGA_Molecule_calc_df
//{
//    struct d2idx_type
//    {
//	d2idx_type() : d2(0), i(0), j(0) { }
//	d2idx_type(int nd2, int ni, int nj) : d2(nd2), i(ni), j(nj) { }
//	int d2, i, j;
//    };
//    bool operator<(const d2idx_type& lhs, const d2idx_type& rhs)
//    {
//	return lhs.d2 < rhs.d2;
//    }
//}
//
//void Molecule::set_mbad_abadMax() const
//{
//    mbad = abad.sum();
//    abadMax = (NAtoms > 0) ? max(abad.max(), (double) NAtoms) : 0.0;
//}

// recalculate everything from scratch
void Molecule::Recalculate()
{
    if (NAtoms() >= ss->NDist)
    {
	cerr << "E: molecule too large" << endl;
	throw InvalidMolecule();
    }
    // destroy all pairs
    typedef list<Pair_t*>::iterator LPPit;
    for (LPPit ii = pairs.begin(); ii != pairs.end(); ++ii)
    {
	delete *ii;
    }
    pairs.clear();
    // molecule parameters
    max_abad = -1;
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
	    pairs.push_back(new Pair_t(this, *ai, *aj));
	}
    }
}

//double Molecule::ABadness(int i) const
//{
//    if (!cached) calc_db();
//    return abad[i];
//}
//
//double Molecule::AFitness(int i) const
//{
//    // this will update abadMax if necessary
//    double badness_i = ABadness(i);
//    return abadMax - badness_i;
//}
//
int Molecule::Badness() const
{
    return badness;
}

int Molecule::Fitness() const
{
    // call to MaxABadness() will update max_abad if necessary
    return NAtoms()*MaxABadness() - Badness();
}

bool comp_Atom_Badness(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.Badness() < rhs.Badness();
}

bool comp_Atom_h(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.h < rhs.h;
}

bool comp_Atom_k(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.Badness() < rhs.Badness();
}

bool comp_Atom_l(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.l < rhs.l;
}

int Molecule::MaxABadness() const
{
    if (max_abad < 0)
    {
	max_abad = (NAtoms() == 0) ?  0 : max( NAtoms()-1, 
		max_element(atoms.begin(), atoms.end(),
		    comp_Atom_Badness)-> Badness() );
    }
    return max_abad;
}

//double Molecule::ABadnessAt(int nh, int nk) const
//{
//    if (NAtoms == ss->NAtoms)
//    {
//	cerr << "E: molecule too large, in Molecule::ABadnessAt()" << endl;
//	throw InvalidMolecule();
//    }
//    if (!cached) calc_db();
//    valarray<int> nd2(NAtoms);
//    for (int dhi, dki, i = 0; i < NAtoms; ++i)
//    {
//	dhi = h[i] - nh;
//	dki = k[i] - nk;
//	nd2[i] = dhi*dhi + dki*dki;
//    }
//    sort(&nd2[0], &nd2[nd2.size()]);
//    double nbad = 0.0;
//    vector<list>::iterator idif = ssdIdxFree.begin();
//    for (int *pnd2 = &(nd2[0]); pnd2 != &(nd2[nd2.size()]); ++pnd2)
//    {
//	for ( ; idif != ssdIdxFree.end() && ss->d2hi[*idif] < *pnd2;
//		++idif ) {};
//	// we found a baddie if we reached the end of ssdIdxFree table
//	// or if the current distance is smaller than d2lo
//	if (idif == ssdIdxFree.end() || *pnd2 < ss->d2lo[*idif])
//	{
//	    nbad++;
//	    if (*pnd2 < ss->d2lo[0])  nbad += NAtoms;
//	}
//	// otherwise it is a matching distance
//	else
//	{
//	    ++idif;
//	}
//    }
//    return nbad;
//}
//
//double Molecule::MBadnessWith(const Molecule& M) const
//{
//    if (NAtoms+M.NAtoms > ss->NAtoms)
//    {
//	cerr << "E: joined molecule too large" << endl;
//	throw InvalidMolecule();
//    }
//    if (!cached) calc_db();
//    if (!M.cached) M.calc_db();
//    valarray<int> nd2(NAtoms*M.NAtoms);
//    for (int dh, dk, id2 = 0, i = 0; i < NAtoms; ++i)
//    {
//	for (int j = 0; j < M.NAtoms; ++j)
//	{
//	    dh = h[i] - M.h[j];
//	    dk = k[i] - M.k[j];
//	    nd2[id2++] = dh*dh + dk*dk;
//	}
//    }
//    sort(&nd2[0], &nd2[nd2.size()]);
//    double mbadwith = MBadness() + M.MBadness();
//    vector<int>::iterator idif = ssdIdxFree.begin();
//    for (int *pnd2 = &(nd2[0]); pnd2 != &(nd2[nd2.size()]); ++pnd2)
//    {
//	for ( ; idif != ssdIdxFree.end() && ss->d2hi[*idif] < *pnd2;
//		++idif ) {};
//	// we found a baddie if we reached the end of ssdIdxFree table
//	// or if the current distance is smaller than d2lo
//	if (idif == ssdIdxFree.end() || *pnd2 < ss->d2lo[*idif])
//	{
//	    mbadwith += 2.0;
//	}
//	// otherwise it is a matching distance
//	else
//	{
//	    ++idif;
//	}
//    }
//    // now add penalty for being outside the SandSphere
//    return mbad;
//}


//////////////////////////////////////////////////////////////////////////
// Molecule operators
//////////////////////////////////////////////////////////////////////////

Molecule& Molecule::Shift(double dh, double dk, double dl)
{
    for (list<Atom_t>::iterator ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	ai->h += (int) round(dh);
	ai->k += (int) round(dk);
	ai->l += (int) round(dl);
    }
    return *this;
}

Molecule& Molecule::Center()
{
    double avg_h = 0.0, avg_k = 0.0, avg_l = 0.0;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	avg_h += ai->h;
	avg_k += ai->k;
	avg_l += ai->l;
    }
    avg_h /= NAtoms();
    avg_k /= NAtoms();
    avg_l /= NAtoms();
    Shift(-avg_h, -avg_k, -avg_l);
    return *this;
}

//Molecule& Molecule::Rotate(double phi, double h0, double k0)
//{
//    // define rotation matrix
//    double Rz[2][2] = {
//	{ cos(phi),	-sin(phi)  },
//	{ sin(phi), 	 cos(phi)  }
//    };
//    double hshift, kshift;
//    for (int i = 0; i < NAtoms; ++i)
//    {
//	// shift coordinates origin to the center of rotation
//	double hshift = h[i] - h0;
//	double kshift = k[i] - k0;
//	// rotate
//	double hrot = Rz[0][0]*hshift + Rz[0][1]*kshift;
//	double krot = Rz[1][0]*hshift + Rz[1][1]*kshift;
//	// shift back
//	h[i] = (int) round(hrot + h0);
//	k[i] = (int) round(krot + k0);
//    }
//    return *this;
//}
//
//Molecule& Molecule::Part(const Molecule& M, const int cidx)
//{
//    h.resize(1, M.h[cidx]);
//    k.resize(1, M.k[cidx]);
//    fix_size();
//    return *this;
//}
//
//Molecule& Molecule::Part(const Molecule& M, const list<int>& cidx)
//{
//    vector<int> *h_new, *k_new;
//    if (&M == this)
//    {
//	h_new = new vector<int>;
//	k_new = new vector<int>;
//    }
//    else
//    {
//	h_new = &h;
//	k_new = &k;
//    }
//    h_new->resize(cidx.size());
//    k_new->resize(cidx.size());
//    int idx = 0;
//    for ( list<int>::const_iterator li = cidx.begin();
//	    li != cidx.end(); ++li )
//    {
//	if (*li < 0 || *li >= M.NAtoms)
//	{
//	    throw range_error("in Molecule::Part()");
//	}
//	(*h_new)[idx] = M.h[*li];
//	(*k_new)[idx] = M.k[*li];
//	++idx;
//    }
//    if (&M == this)
//    {
//	h = *h_new;
//	k = *k_new;
//	delete h_new, k_new;
//    }
//    fix_size();
//    return *this;
//}

Molecule& Molecule::Pop(list<Atom_t>::iterator ai)
{
    // delete all pairs that refer to *ai
    Atom_t *ap = &(*ai);
    typedef list<Pair_t*>::iterator LPPit;
    for (LPPit ii = pairs.begin(); ii != pairs.end(); ++ii)
    {
	Pair_t* pp = *ii;
	if (pp->atom1 == ap || pp->atom2 == ap)
	{
	    delete *ii;
	    pairs.erase(ii);
	}
    }
    atoms.erase(ai);
    // uncache max atom badness
    max_abad = -1;
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

//Molecule& Molecule::Pop(const Molecule& M, const int cidx)
//{
//    if (cidx < 0 || cidx >= M.NAtoms)
//    {
//	throw range_error("in Molecule::Pop(list<int>)");
//    }
//    if (this != &M)  *this = M;
//    h.erase(h.begin() + cidx);
//    k.erase(k.begin() + cidx);
//    fix_size();
//    return *this;
//}
//
//Molecule& Molecule::Pop(const Molecule& M, const list<int>& cidx)
//{
//    if (cidx.size() == 0)
//    {
//	if (this != &M)  *this = M;
//	return *this;
//    }
//    list<int> sidx(cidx);
//    sidx.sort();
//    if (sidx.front() < 0 || sidx.back() >= M.NAtoms)
//    {
//	throw range_error("in Molecule::Pop(list<int>)");
//    }
//    sidx.push_back(M.NAtoms);
//    vector<int> h_new, k_new;
//    int j = 0;
//    for ( list<int>::iterator li = sidx.begin();
//	    li != sidx.end(); ++li )
//    {
//	for (; j < *li; ++j)
//	{
//	    h_new.push_back(M.h[j]);
//	    k_new.push_back(M.k[j]);
//	}
//	j = *li + 1;
//    }
//    h = h_new;
//    k = k_new;
//    fix_size();
//    return *this;
//}

Molecule& Molecule::Clear()
{
    // pairs must be destroyed before atoms;
    typedef list<Pair_t*>::iterator LPPit;
    for (LPPit ii = pairs.begin(); ii != pairs.end(); ++ii)
    {
	delete *ii;
    }
    pairs.clear();
    atoms.clear();
    // uncache max atom badness
    max_abad = -1;
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

Molecule& Molecule::Add(int nh, int nk, int nl)
{
    Add(Atom_t(nh, nk, nl));
    return *this;
}

Molecule& Molecule::Add(Atom_t atom)
{
    if (NAtoms() == max_NAtoms)
    {
	cerr << "E: molecule too large" << endl;
	throw InvalidMolecule();
    }
    list<Atom_t>::iterator this_atom, ai;
    this_atom = atoms.insert(atoms.end(), atom);
    this_atom->ResetBadness();
    // uncache max atom badness
    max_abad = -1;
    for (ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	if (ai == this_atom)  continue;
	pairs.push_back(new Pair_t(this, *ai, *this_atom));
    }
    return *this;
}

//Molecule& Molecule::Add(Molecule& M)
//{
//    for (int i = 0; i < M.NAtoms; ++i)
//    {
//	h.push_back(M.h[i]);
//	k.push_back(M.k[i]);
//    }
//    fix_size();
//    return *this;
//}
//
//Molecule& Molecule::Add(int nh, int nk)
//{
//    h.push_back(nh);
//    k.push_back(nk);
//    fix_size();
//    return *this;
//}
//
//Molecule& Molecule::MoveAtomTo(int idx, int nh, int nk)
//{
//    if (idx >= NAtoms)
//    {
//	throw range_error("in Molecule::MoveAtomTo()");
//    }
//    UnCache();
//    h[idx] = nh;
//    k[idx] = nk;
//    return *this;
//}
//
//struct Molecule::badness_at
//{
//    badness_at() : h(0), k(0)
//    {
//	numeric_limits<double> double_info;
//	abad = double_info.max();
//    }
//    badness_at(int nh, int nk, double nbad) :
//	h(nh), k(nk), abad(nbad) { }
//    int h, k;
//    double abad;
//};
//bool operator<(
//	const Molecule::badness_at& lhs, const Molecule::badness_at& rhs
//	)
//{
//    return lhs.abad < rhs.abad;
//}
//
//list<Molecule::badness_at> Molecule::find_good_distances(
//	int trials, const vector<int>& ssd_idx
//	)
//{
//    if (NAtoms > ss->NDist)
//    {
//	cerr << "E: molecule too large for finding a new position" << endl;
//	throw InvalidMolecule();
//    }
//    // prepare a return list:
//    list<badness_at> retlist;
//    // prepare discrete RNG
//    double afit[NAtoms];
//    for (int i = 0; i != NAtoms; ++i)
//    {
//	afit[i] = AFitness(i);
//    }
//    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(NAtoms, afit);
//    for (int i = 0; i < trials; ++i)
//    {
//	// pick one atom and free distance
//	int a1 = gsl_ran_discrete(BGA::rng, table);
//	int idx = gsl_rng_uniform_int(BGA::rng, ssd_idx.size());
//	double radius = ss->d[ssd_idx[idx]];
//	double phi = 2.0*M_PI*gsl_rng_uniform(BGA::rng);
//	int nh = h[a1] + (int)round(radius * cos(phi));
//	int nk = k[a1] + (int)round(radius * sin(phi));
//	retlist.push_back( badness_at(nh, nk, ABadnessAt(nh, nk)) );
//    }
//    gsl_ran_discrete_free(table);
//    return retlist;
//}
//
//list<Molecule::badness_at> Molecule::find_good_triangles(
//	int trials, const vector<int>& ssd_idx
//	)
//{
//    // try to generate 2 triangles with good distances, this may not always
//    // work, so this returns a list, which is either empty, or has even
//    // number of entries
//    if (NAtoms > ss->NDist)
//    {
//	cerr << "E: molecule too large for finding a new position" << endl;
//	throw InvalidMolecule();
//    }
//    else if (NAtoms < 2)
//    {
//	cerr << "E: molecule too small, triangulation not possible" << endl;
//	throw InvalidMolecule();
//    }
//    // prepare a return list:
//    list<badness_at> retlist;
//    // prepare discrete RNG
//    double afit[NAtoms];
//    for (int i = 0; i != NAtoms; ++i)
//    {
//	afit[i] = AFitness(i);
//    }
//    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(NAtoms, afit);
//    for (int nt = 0; nt < trials; ++nt)
//    {
//	// pick first atom and free distance
//	int a1 = gsl_ran_discrete(BGA::rng, table);
//	int a2 = gsl_ran_discrete(BGA::rng, table);
//	for (int i = 0; i < 5 && a1 == a2; ++i)
//	{
//	    a2 = gsl_ran_discrete(BGA::rng, table);
//	}
//	if (a1 == a2)
//	{
//	    a2 = (a1 + 1) % NAtoms;
//	}
//	int idf1 = gsl_rng_uniform_int(BGA::rng, ssd_idx.size());
//	int idf2 = gsl_rng_uniform_int(BGA::rng, ssd_idx.size()-1) + 1;
//	idf2 = (idf2 + idf1) % ssd_idx.size();
//	if (idf1 == idf2) throw(runtime_error("idf1 == idf2"));
//	double r13 = ss->d[ssd_idx[idf1]];
//	double r23 = ss->d[ssd_idx[idf2]];
//	double r12 = dist(a1, a2);
//	// is triangle [a1 a2 a3] possible?
//	if (r12 < 1.0 || r13 + r23 < r12 || fabs(r13 - r23) > r12) continue;
//	// here we can construct 2 triangles
//	double longdir[2] = {  (h[a2]-h[a1])/r12, (k[a2]-k[a1])/r12 };
//	double perpdir[2] = { -longdir[1], 	  longdir[0] };
//	double xlong = (r13*r13 + r12*r12 - r23*r23) / (2.0*r12);
//	double xperp = sqrt(r13*r13 - xlong*xlong);
//	int nh1 = h[a1] + (int) round(xlong*longdir[0] + xperp*perpdir[0]);
//	int nk1 = k[a1] + (int) round(xlong*longdir[1] + xperp*perpdir[1]);
//	// store the result
//	badness_at res1(nh1, nk1, ABadnessAt(nh1, nk1));
//	retlist.push_back(res1);
//	// jump out if it is a good location or if we did enough trials
//	/* debug:
//	   cout << nt << " a1 = " << a1 << " a2 = " << a2 << endl;
//	   cout << nt << " res1.abad = " << res1.abad << endl;
//	   cout << nt << " res1.h = " << res1.h <<
//	       " res1.k = " << res1.k << endl;
//	   */
//	if (res1.abad == 0.0 || ++nt == trials)  break;
//	// 2nd triangle:
//	int nh2 = h[a1] + (int) round(xlong*longdir[0] - xperp*perpdir[0]);
//	int nk2 = k[a1] + (int) round(xlong*longdir[1] - xperp*perpdir[1]);
//	// store the result
//	badness_at res2(nh2, nk2, ABadnessAt(nh2, nk2));
//	retlist.push_back(res2);
//	/* debug:
//	   cout << nt << " res2.abad = " << res2.abad << endl;
//	   cout << nt << " res2.h = " << res2.h <<
//	       " res2.k = " << res2.k << endl;
//	   cout << nt <<
//	   " longdir=[" << longdir[0] << ' ' << longdir[1] << "]," <<
//	   " perpdir=[" << perpdir[0] << ' ' << perpdir[1] << "]," <<
//	   " xlong=" << xlong << "," <<
//	   " xperp=" << xperp << "," <<
//	   " r12=" << r12 << "," <<
//	   " r13=" << r13 << "," <<
//	   " r23=" << r23 << "," <<
//	   endl;
//	   */
//	if (res2.abad == 0.0 || ++nt == trials)  break;
//    }
//    gsl_ran_discrete_free(table);
//    return retlist;
//}
//
//list<Molecule::badness_at> Molecule::find_good_triangles2(
//	int trials, const vector<int>& ssd_idx
//	)
//{
//    // try to generate 2 triangles with good distances, this may not always
//    // work, so this returns a list, which is either empty, or has even
//    // number of entries
//    if (NAtoms > ss->NDist)
//    {
//	cerr << "E: molecule too large for finding a new position" << endl;
//	throw InvalidMolecule();
//    }
//    else if (NAtoms < 2)
//    {
//	cerr << "E: molecule too small, triangulation not possible" << endl;
//	throw InvalidMolecule();
//    }
//    // prepare a return list:
//    list<badness_at> retlist;
//    for (int nt = 0; nt < trials; ++nt)
//    {
//	// choose 2 free indices
//	int idf1 = gsl_rng_uniform_int(BGA::rng, ssd_idx.size());
//	int idf2 = gsl_rng_uniform_int(BGA::rng, ssd_idx.size()-1) + 1;
//	idf2 = (idf2 + idf1) % ssd_idx.size();
//	if (idf1 == idf2) throw(runtime_error("idf1 == idf2"));
//	double r13 = ss->d[ssd_idx[idf1]];
//	double r23 = ss->d[ssd_idx[idf2]];
//	// loop over all atom pairs
//	for (int a1 = 0; a1 < NAtoms; ++a1)
//	{
//	    for (int a2 = a1+1; a2 < NAtoms; ++a2)
//	    {
//		double r12 = dist(a1, a2);
//		// is triangle [a1 a2 a3] possible?
//		if (r12 < 1.0 || r13 + r23+0.5 < r12 || fabs(r13 - r23) > r12+0.5)
//		    continue;
//		// here we can construct 2 triangles
//		double longdir[2] = {  (h[a2]-h[a1])/r12, (k[a2]-k[a1])/r12 };
//		double perpdir[2] = { -longdir[1], 	  longdir[0] };
//		double xlong = (r13*r13 + r12*r12 - r23*r23) / (2.0*r12);
//		double xperp2 = r13*r13 - xlong*xlong;
//		double xperp = (xperp2 > 0.0) ? sqrt(r13*r13 - xlong*xlong) : 0;
//		int nh1 = h[a1] + (int) round(xlong*longdir[0] + xperp*perpdir[0]);
//		int nk1 = k[a1] + (int) round(xlong*longdir[1] + xperp*perpdir[1]);
//		// store the result
//		badness_at res1(nh1, nk1, abad[a1]+abad[a2]+ABadnessAt(nh1, nk1));
//		retlist.push_back(res1);
//		// jump out if it is a good location or if we did enough trials
//		/* debug:
//		   cout << nt << " a1 = " << a1 << " a2 = " << a2 << endl;
//		   cout << nt << " res1.abad = " << res1.abad << endl;
//		   cout << nt << " res1.h = " << res1.h <<
//		   " res1.k = " << res1.k << endl;
//		   */
//		// 2nd triangle:
//		int nh2 = h[a1] + (int) round(xlong*longdir[0] - xperp*perpdir[0]);
//		int nk2 = k[a1] + (int) round(xlong*longdir[1] - xperp*perpdir[1]);
//		// store the result
//		badness_at res2(nh2, nk2, abad[a1]+abad[a2]+ABadnessAt(nh2, nk2));
//		retlist.push_back(res2);
//		/* debug:
//		   cout << nt << " res2.abad = " << res2.abad << endl;
//		   cout << nt << " res2.h = " << res2.h <<
//		   " res2.k = " << res2.k << endl;
//		   cout << nt <<
//		   " longdir=[" << longdir[0] << ' ' << longdir[1] << "]," <<
//		   " perpdir=[" << perpdir[0] << ' ' << perpdir[1] << "]," <<
//		   " xlong=" << xlong << "," <<
//		   " xperp=" << xperp << "," <<
//		   " r12=" << r12 << "," <<
//		   " r13=" << r13 << "," <<
//		   " r23=" << r23 << "," <<
//		   endl;
//		   */
//	    }
//	}
//    }
//    return retlist;
//}
//
//Molecule& Molecule::Evolve(int trials)
//{
//    if (NAtoms == ss->NAtoms)
//    {
//	cerr << "E: full-sized molecule cannot Evolve()" << endl;
//	throw InvalidMolecule();
//    }
//    // evolution is trivial for empty or 1-atom molecule
//    switch (NAtoms)
//    {
//	case 0:
//	    Add(0, 0);
//	    return *this;
//	case 1:
//	    // make sure ssdIdxFree is updated
//	    if (!cached) calc_db();
//	    list<badness_at> ba = find_good_distances(1, ssdIdxFree);
//	    Add(ba.front().h, ba.front().k);
//	    Center();
//	    return *this;
//    }
//    // make sure ssdIdxFree is updated
//    if (!cached) calc_db();
//    // here we can be sure that NAtoms >= 2
//    list<badness_at> trials_log = find_good_distances(trials/3, ssdIdxFree);
//    // if we did not find a good place, let's try triangulation
//    if (trials_log.back().abad != 0.0)
//    {
//	list<badness_at> fgt =
//	    find_good_triangles2(1, ssdIdxFree);
////	    find_good_triangles(trials-trials_log.size(), ssdIdxFree);
//	trials_log.insert(trials_log.end(), fgt.begin(), fgt.end());
//    }
//    // purge trials_log to minimal elements
//    badness_at best = *min_element(trials_log.begin(), trials_log.end());
//    trials_log.remove_if( bind1st(less<badness_at>(),best) );
//    int idx = gsl_rng_uniform_int(BGA::rng, trials_log.size());
//    list<badness_at>::iterator lbi = trials_log.begin(); advance(lbi, idx);
//    best = *lbi;
//    Add(best.h, best.k);
//    Center();
//    return *this;
//}
//
//namespace BGA_Molecule_MateWith
//{
//    list<int> random_wt_choose(int N, const double *p, int Np)
//    {
//	list<int> retlst;
//	if (N > Np)
//	{
//	    throw(range_error("too many items to choose"));
//	}
//	// check trivial cases
//	else if (N == 0)
//	{
//	    return retlst;
//	}
//	else if (N == Np)
//	{
//	    for (int i = 0; i < Np; ++i)
//	    {
//		retlst.push_back(i);
//	    }
//	    return retlst;
//	}
//	if ( p+Np != find_if(p, p+Np, bind2nd(less<double>(),0.0)) )
//	{
//	    throw(runtime_error("negative choice probability"));
//	}
//	// now we need to do some real work
//	double prob[Np];
//	copy(p, p+Np, prob);
//	// integer encoding
//	int val[Np];
//	for (int i = 0; i != Np; ++i)  val[i] = i;
//	// cumulative probability
//	double cumprob[Np];
//	int Nprob = Np;
//	// main loop
//	for (int i = 0, Nprob = Np; i != N; ++i, --Nprob)
//	{
//	    // calculate cumulative probability
//	    partial_sum(prob, prob+Nprob, cumprob);
//	    // if all probabilities are 0.0, set them to equal value
//	    if (cumprob[Nprob-1] == 0.0)
//	    {
//		for (int j = 0; j != Nprob; ++j)
//		{
//		    prob[j] = 1.0;
//		    cumprob[j] = (j+1.0)/Nprob;
//		}
//	    }
//	    // otherwise we can normalize cumprob
//	    else
//	    {
//		for (double *pcp = cumprob; pcp != cumprob+Nprob; ++pcp)
//		{
//		    *pcp /= cumprob[Nprob-1];
//		}
//	    }
//	    // now let's do binary search on cumprob:
//	    double r = gsl_rng_uniform(BGA::rng);
//	    double *pcp = upper_bound(cumprob, cumprob+Nprob, r);
//	    int idx = pcp - cumprob;
//	    retlst.push_back(val[idx]);
//	    // overwrite this element with the last number
//	    prob[idx] = prob[Nprob-1];
//	    val[idx] = val[Nprob-1];
//	}
//	return retlst;
//    }
//
//    bool compare_MFitness(Molecule& lhs, Molecule& rhs)
//    {
//	return lhs.MFitness() < rhs.MFitness();
//    }
//}
//
//Molecule& Molecule::Degenerate(int Npop)
//{
//    using namespace BGA_Molecule_MateWith;
//    Npop = min(NAtoms, Npop);
//    if (Npop == 0)  return *this;
//    // make sure valarray abad is updated
//    if (!cached) calc_db();
//    list<int> ipop = random_wt_choose(Npop, &(abad[0]), NAtoms);
//    Pop(ipop);
//    Center();
//    return *this;
//}
//
//Molecule Molecule::MateWith(const Molecule& Male, int trials)
//{
//    using namespace BGA_Molecule_MateWith;
//    if (NAtoms != Male.NAtoms)
//    {
//	cerr << "E: mating of unequal molecules" << endl;
//	throw InvalidMolecule();
//    }
//    double patom_f[NAtoms];
//    for (int i = 0; i < NAtoms; ++i)  patom_f[i] = AFitness(i);
//    double patom_m[Male.NAtoms];
//    for (int i = 0; i < Male.NAtoms; ++i)  patom_m[i] = Male.AFitness(i);
//    list<Molecule> children;
//    for (int i = 0; i < trials; ++i)
//    {
//	Molecule egg(ss);
//	list<int> egg_atoms =
//	    random_wt_choose(NAtoms, patom_f, NAtoms/2);
//	egg.Part(*this, egg_atoms);
//	Molecule sperm(ss);
//	list<int> sperm_atoms =
//	    random_wt_choose(Male.NAtoms, patom_m, NAtoms - egg.NAtoms);
//	sperm.Part(Male, sperm_atoms);
//	egg.mount(sperm);
//	children.push_back(egg);
//	if (egg.MBadness() == 0.0)  break;
//    }
//    // find the best child
//    Molecule *best;
//    // maybe the last child is a prodigy
//    if (children.back().MBadness() == 0.0)
//    {
//	best = &children.back();
//    }
//    else
//    {
//	best = &( *min_element(children.begin(), children.end(),
//		compare_MFitness) );
//    }
//    best->Center();
//    return *best;
//}
//
//Molecule& Molecule::mount(Molecule& sperm)
//{
//    using namespace BGA_Molecule_MateWith;
//    if (NAtoms < sperm.NAtoms)
//    {
//	cerr << "E: sperm molecule is larger than egg" << endl;
//	throw InvalidMolecule();
//    }
//    // mounting is trivial for empty or 1-atom molecules
//    if (NAtoms == 0 || sperm.NAtoms ==0)
//    {
//	return *this;
//    }
//    else if (NAtoms == 1)
//    {
//	// here sperm.NAtoms must be equal 1
//	// make sure ssdIdxFree is updated
//	if (!cached) calc_db();
//	list<badness_at> ba = find_good_distances(1, ssdIdxFree);
//	Add(ba.front().h, ba.front().k);
//	Center();
//	return *this;
//    }
//    else if (NAtoms == 2 && sperm.NAtoms == 1)
//    {
//	Evolve();
//	return *this;
//    }
//    // here we can be sure that NAtoms >= 2, sperm.NAtoms >= 2
//    // make sure ssdIdxFree is updated for egg and sperm
//    if (!cached) calc_db();
//    if (!sperm.cached) sperm.calc_db();
//    // calculate probabilities of choosing an egg free distance,
//    // free distances common to egg and sperm have 5 times higher probability
//    const double common_prob = 5.0;
//    double pidx[ssdIdxFree.size()];
//    double *pp = pidx;
//    for (vector<int>::iterator
//	    eidx = ssdIdxFree.begin(), sidx = sperm.ssdIdxFree.begin();
//	    eidx != ssdIdxFree.end(); ++eidx, ++pp)
//    {
//	for (; sidx != sperm.ssdIdxFree.end() && *sidx < *eidx; ++sidx) { }
//	// no more common distances if sidx was used
//	if (sidx != sperm.ssdIdxFree.end() && *sidx == *eidx)
//	{
//	    *pp = common_prob;
//	}
//	else
//	{
//	    *pp = 1.0;
//	}
//    }
//    // let's start with random mount:
//    int idf = random_wt_choose(ssdIdxFree.size(), pidx, 1).front();
//    double e_afit[NAtoms];
//    for (int i = 0; i < NAtoms; ++i)  e_afit[i] = AFitness(i);
//    double s_afit[sperm.NAtoms];
//    for (int i = 0; i < sperm.NAtoms; ++i)  s_afit[i] = sperm.AFitness(i);
//    int ea =  random_wt_choose(NAtoms, e_afit, 1).front();
//    int sa =  random_wt_choose(sperm.NAtoms, s_afit, 1).front();
//    double omega = 2.0*M_PI*gsl_rng_uniform(BGA::rng);
//    Molecule sp1(sperm);
//    sp1.Rotate(omega, sp1.h[sa], sp1.k[sa]);
//    double radius = ss->d[ssdIdxFree[idf]];
//    double phi = 2.0*M_PI*gsl_rng_uniform(BGA::rng);
//    double shn = h[ea] + radius*cos(phi);
//    double skn = k[ea] + radius*sin(phi);
//    sp1.Shift(shn-sp1.h[sa], skn-sp1.k[sa]);
//    double ch1_badness = MBadnessWith(sp1);
//    // save child molecule
//    Molecule ch1(*this); ch1.Add(sp1);
//    // next try a sophisticated triangular mount:
//    list<int> lidf;
//    lidf = random_wt_choose(NAtoms, e_afit, 2);
//
//}


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
	read_token("NAtoms", NAtoms) &&
	read_token("delta", delta);
    if (!state)
    {
	return;
    }
    if (fmt == "grid")
	format = GRID;
    else if (fmt == "xy")
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

istream& Molecule::ReadGrid(istream& fid)
{
    // read values to integer vector vhkl
    string header;
    vector<int> vhkl;
    bool result = read_header(fid, header) && read_data(fid, vhkl);
    if (!result) return fid;
    // parse header
    double vhkl_scale = 1.0;
    int vhkl_NAtoms = vhkl.size()/3;
    ParseHeader ph(header);
    if (ph)
    {
	vhkl_scale = ph.delta / ss->delta;
	if ( vhkl_NAtoms != ph.NAtoms )
	{
	    cerr << "E: " << opened_file << ": expected " << ph.NAtoms <<
		" atoms, read " << vhkl_NAtoms << endl;
	    throw IOError();
	}
    }
    // check if all coordinates have been read
    if ( vhkl.size() % 3 )
    {
	cerr << "E: " << opened_file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    for (int i = 0; i < vhkl.size(); i += 3)
    {
	int h = vhkl[i+0];
	int k = vhkl[i+1];
	int l = vhkl[i+2];
	Add(Atom_t(h, k, l));
    }
    return fid;
}

bool Molecule::ReadGrid(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    opened_file = file;
    bool result = ReadGrid(fid);
    opened_file.clear();
    fid.close();
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
	int h = (int) round(vxyz[i+0] / ss->delta);
	int k = (int) round(vxyz[i+1] / ss->delta);
	int l = (int) round(vxyz[i+2] / ss->delta);
	Add(Atom_t(h, k, l));
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
    // open file for writing
    ofstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to write to '" << file << "'" << endl;
	throw IOError();
    }
    bool result = (fid << M);
    fid.close();
    return result;
}

bool Molecule::WriteGrid(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtGrid();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
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

Molecule& Molecule::OutFmtGrid()
{
    output_format = GRID;
    return *this;
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
	case M.GRID:
	    result = M.ReadGrid(fid);
	    break;
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
	case M.GRID:
	    fid << "# BGA molecule format = grid" << endl;
	    fid << "# NAtoms = " << M.NAtoms() << endl;
	    fid << "# delta = " << M.ss->delta << endl;
	    for (LAit ai = afirst; ai != alast; ++ai)
	    {
		fid << ai->h << '\t' << ai->k << '\t' << ai->l << endl;
	    }
	    break;
	case M.XYZ:
	    fid << "# BGA molecule format = xy" << endl;
	    fid << "# NAtoms = " << M.NAtoms() << endl;
	    fid << "# delta = " << M.ss->delta << endl;
	    for (LAit ai = afirst; ai != alast; ++ai)
	    {
		fid <<
		    M.ss->delta * ai->h << '\t' <<
		    M.ss->delta * ai->k << '\t' <<
		    M.ss->delta * ai->l << endl;
	    }
	    break;
	case M.ATOMEYE:
	    double xyz_lo = 0.0;
	    double xyz_hi = 1.0;
	    double xyz_range = xyz_hi - xyz_lo;
	    if (M.NAtoms() > 0)
	    {
		const double scale = 1.01*M.ss->delta;
		double xyz_extremes[8] = {
		    -M.ss->dmax,
		    scale * min_element(afirst, alast, comp_Atom_h)->h,
		    scale * min_element(afirst, alast, comp_Atom_k)->k,
		    scale * min_element(afirst, alast, comp_Atom_l)->l,
		    M.ss->dmax,
		    scale * max_element(afirst, alast, comp_Atom_h)->h,
		    scale * max_element(afirst, alast, comp_Atom_k)->k,
		    scale * max_element(afirst, alast, comp_Atom_l)->l,
		};
		xyz_lo = *min_element(xyz_extremes, xyz_extremes+8);
		xyz_hi = *max_element(xyz_extremes, xyz_extremes+8);
		xyz_range = xyz_hi - xyz_lo;
	    }
	    fid << "# BGA molecule format = atomeye" << endl;
	    fid << "# NAtoms = " << M.NAtoms() << endl;
	    fid << "# delta = " << M.ss->delta << endl;
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
		    (ai->h * M.ss->delta - xyz_lo) / xyz_range << " " <<
		    (ai->k * M.ss->delta - xyz_lo) / xyz_range << " " <<
		    (ai->l * M.ss->delta - xyz_lo) / xyz_range << " " <<
		    ai->Badness() << endl;
	    }
	    break;
    }
    return fid;
}

void Molecule::PrintBadness()
{
    // call to MBadness() will update abadMax if necessary
    cout << "MBadness() = " << Badness() << endl;
    cout << "ABadness() =";
    double mab = MaxABadness();
    bool marked = false;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	cout << ' ' << ai->Badness();
	if (!marked && ai->Badness() == mab)
	{
	    cout << '*';
	    marked = true;
	}
    }
    cout << endl;
}

void Molecule::PrintFitness()
{
    cout << "MFitness() = " << Fitness() << endl;
    cout << "AFitness() =";
    double mab = MaxABadness();
    bool marked = false;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	cout << ' ' << mab - ai->Badness();
	if (!marked && ai->Badness() == mab)
	{
	    cout << '*';
	    marked = true;
	}
    }
    cout << endl;
}

//////////////////////////////////////////////////////////////////////////
//// Population definitions
//////////////////////////////////////////////////////////////////////////
//
//void Population::init()
//{
//    // check whether all members use the same SandSphere
//    /*
//    iterator ibad =
//	adjacent_find(begin(), end(), molecule_SandSpheres_differ);
//    if (ibad != end())
//    {
//	cerr << "E: population contains molecules with different sandsphere" <<
//	    endl;
//	throw InvalidPopulation();
//    }
//    */
//}
//
//vector<Couple> Population::FindCouples(int NCouples)
//{
//    double fitness[size()];
//    double *pf = fitness;
//    for (iterator i = begin(); i != end(); ++i)
//    {
//	*(pf++) = i->MFitness();
//    }
//    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(size(), fitness);
//    vector<Couple> v(NCouples);
//    typedef vector<Couple>::iterator Cit;
//    for (Cit i = v.begin(); i != v.end(); ++i)
//    {
//	i->Male = &at( gsl_ran_discrete(BGA::rng, table) );
//	i->Female = &at( gsl_ran_discrete(BGA::rng, table) );
//    }
//    gsl_ran_discrete_free(table);
//    return v;
//}
