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
	(*i)->UnCache();
    }
}

double SandSphere::GridTol()
{
    return vGridTol;
}


////////////////////////////////////////////////////////////////////////
// Molecule definitions
////////////////////////////////////////////////////////////////////////

Molecule::Molecule(SandSphere *SS) : ss(SS)
{
    init();
}

Molecule::Molecule(SandSphere *SS,
	int s, int *ph, int *pk) : ss(SS)
{
    h.resize(s);
    k.resize(s);
    for (int i = 0; i < s; ++i)
    {
	h[i] = ph[i];
	k[i] = pk[i];
    }
    init();
}

Molecule::Molecule(SandSphere *SS,
	const vector<int>& vh, const vector<int>& vk) : ss(SS)
{
    h = vh;
    k = vk;
    init();
}

Molecule::Molecule(SandSphere *SS,
	int s, double *px, double *py) : ss(SS)
{
    h.resize(s);
    k.resize(s);
    for (int i = 0; i < s; ++i)
    {
	h[i] = (int) round(px[i] / ss->delta);
	k[i] = (int) round(py[i] / ss->delta);
    }
    init();
}

Molecule::Molecule(SandSphere *SS,
	const vector<double>& vx, const vector<double>& vy) : ss(SS)
{
    h.resize(vx.size());
    k.resize(vx.size());
    for (size_t i = 0; i < h.size(); ++i)
    {
	h[i] = (int) round(vx[i] / ss->delta);
	k[i] = (int) round(vy[i] / ss->delta);
    }
    init();
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
    // data storage
    if (ss != M.ss)
    {
	ss->molecules.remove(this);
	ss = M.ss;
	ss->molecules.push_back(this);
    }
    // h, k assignment must preceed fix_size()
    h = M.h;			// x-coordinates
    k = M.k;			// y-coordinates
    // parameters
    if (NAtoms != M.NAtoms)
    {
	// this sets NAtoms, NDist, resizes all valarray's
	fix_size();
    }
    // badness evaluation
    cached = M.cached;
    if (M.cached)
    {
	abad = M.abad;		// individual atom badnesses
	abadMax = M.abadMax;	// maximum atom badness
	mbad = M.mbad;		// molecular badness
	d2 = M.d2;		// sorted table of squared distances
	ssdIdxFree = M.ssdIdxFree;	// available elements in ss.dist
    }
    // IO helpers
    output_format = M.output_format;
    opened_file = M.opened_file;
    return *this;
}

void Molecule::init()
{
    ss->molecules.push_back(this);
    MaxAtoms = ss->NAtoms;
    // check coordinate sizes
    OutFmtGrid();
    UnCache();
    fix_size();
}

void Molecule::fix_size()
{
    UnCache();
    if (h.size() != k.size())
    {
	cerr << "E: invalide coordinate vectors" << endl;
	throw InvalidMolecule();
    }
    NAtoms = h.size();
    NDist  = NAtoms*(NAtoms-1)/2;
    if (NAtoms > ss->NAtoms)
    {
	cerr << "E: molecule too large" << endl;
	throw InvalidMolecule();
    }
    abad.resize(NAtoms, 0.0);
    abadMax = 0.0;
    d2.resize(NDist, 0);
//    ssdIdxUsed.clear();
    ssdIdxFree.clear();
}

Molecule::~Molecule()
{
    // debug: cout << "ss->molecules.size() = " << ss->molecules.size() << endl;
    ss->molecules.remove(this);
}


////////////////////////////////////////////////////////////////////////
// Molecule badness/fitness evaluation
////////////////////////////////////////////////////////////////////////

double Molecule::dist(const int& i, const int& j) const
{
    return sqrt(1.0*dist2(i, j));
}

double Molecule::out_penalty(int nh, int nk) const
{
    double Rout = sqrt(nh*nh + nk*nk + 0.0) - ss->gridmax;
    return (Rout > 0.0) ? Rout : 0.0;
}

namespace BGA_Molecule_calc_df
{
    struct d2idx_type
    {
	d2idx_type() : d2(0), i(0), j(0) { }
	d2idx_type(int nd2, int ni, int nj) : d2(nd2), i(ni), j(nj) { }
	int d2, i, j;
    };
    bool operator<(const d2idx_type& lhs, const d2idx_type& rhs)
    {
	return lhs.d2 < rhs.d2;
    }
}

void Molecule::subtract_out_penalty() const
{
    for (int i = 0; cached && i < NAtoms; ++i)
	abad[i] -= out_penalty(h[i], k[i]);
}

void Molecule::add_out_penalty() const
{
    for (int i = 0; cached && i < NAtoms; ++i)
	abad[i] += out_penalty(h[i], k[i]);
}

void Molecule::set_mbad_abadMax() const
{
    mbad = abad.sum();
    abadMax = (NAtoms > 0) ? max(abad.max(), (double) NAtoms) : 0.0;
}

// calculate all distances and atom/molecule badness 
void Molecule::calc_db() const
{
    using namespace BGA_Molecule_calc_df;
    cached = true;
//    ssdIdxUsed.clear();
    ssdIdxFree.clear();
    d2idx_type d2idx[NDist];
    // check if molecule is not too large
    if (NDist > ss->NDist)
    {
	cerr << "E: molecule too large" << endl;
	throw InvalidMolecule();
    }
    // calculate and store distances
    int ij = 0;
    for (int i = 0; i < NAtoms; ++i)
    {
	for (int j = i + 1; j < NAtoms; ++j, ++ij)
	{
	    d2idx[ij].d2 = dist2(i, j);
	    d2idx[ij].i = i;
	    d2idx[ij].j = j;
	}
    }
    sort(d2idx, d2idx+NDist);
    for (ij = 0; ij < NDist; ++ij)
    {
	d2[ij] = d2idx[ij].d2;
    }
    // evaluate abad[i]
    abad = 0.0;
    int d_idx = 0;
    for (ij = 0; ij < NDist; ++ij)
    {
	for(; d_idx < ss->NDist && ss->d2hi[d_idx] < d2idx[ij].d2; ++d_idx)
	{
	    ssdIdxFree.push_back(d_idx);
	}
	// we found a baddie if we reached the end of ss->d2 table
	// or if the current distance is smaller than d2lo
	if (!(d_idx < ss->NDist) || d2idx[ij].d2 < ss->d2lo[d_idx])
	{
	    abad[d2idx[ij].i]++;
	    abad[d2idx[ij].j]++;
	}
	// otherwise it is a matching distance
	else
	{
	    d_idx++;
	}
    }
    // all unused distances from the table are free:
    for(;  d_idx < ss->NDist; ++d_idx)
    {
	ssdIdxFree.push_back(d_idx);
    }
    // now add penalty for being outside the SandSphere
    // this works only when cached == true
    add_out_penalty();
    set_mbad_abadMax();
}

double Molecule::ABadness(int i) const
{
    if (!cached) calc_db();
    return abad[i];
}

double Molecule::AFitness(int i) const
{
    // this will update abadMax if necessary
    double badness_i = ABadness(i);
    return abadMax - badness_i;
}

double Molecule::MBadness() const
{
    if (!cached) calc_db();
    return mbad;
}

double Molecule::MFitness() const
{
    // this will update abadMax if necessary
    double mbadness = MBadness();
    return NAtoms*abadMax - mbadness;
}

double Molecule::ABadnessAt(int nh, int nk) const
{
    if (NAtoms == ss->NAtoms)
    {
	cerr << "E: molecule too large, in Molecule::ABadnessAt()" << endl;
	throw InvalidMolecule();
    }
    if (!cached) calc_db();
    valarray<int> nd2(NAtoms);
    for (int dhi, dki, i = 0; i < NAtoms; ++i)
    {
	dhi = h[i] - nh;
	dki = k[i] - nk;
	nd2[i] = dhi*dhi + dki*dki;
    }
    sort(&nd2[0], &nd2[nd2.size()]);
    double nbad = 0.0;
    list<int>::iterator idif = ssdIdxFree.begin();
    for (int *pnd2 = &(nd2[0]); pnd2 != &(nd2[nd2.size()]); ++pnd2)
    {
	for ( ; idif != ssdIdxFree.end() && ss->d2hi[*idif] < *pnd2;
		++idif ) {};
	// we found a baddie if we reached the end of ssdIdxFree table
	// or if the current distance is smaller than d2lo
	if (idif == ssdIdxFree.end() || *pnd2 < ss->d2lo[*idif])
	{
	    nbad++;
	}
	// otherwise it is a matching distance
	else
	{
	    ++idif;
	}
    }
    // now add penalty for being outside the SandSphere
    nbad += out_penalty(nh, nk);
    return nbad;
}

double Molecule::MBadnessWith(const Molecule& M) const
{
    if (NAtoms+M.NAtoms > ss->NAtoms)
    {
	cerr << "E: joined molecule too large" << endl;
	throw InvalidMolecule();
    }
    if (!cached) calc_db();
    if (!M.cached) M.calc_db();
    valarray<int> nd2(NAtoms*M.NAtoms);
    for (int dh, dk, id2 = 0, i = 0; i < NAtoms; ++i)
    {
	for (int j = 0; j < M.NAtoms; ++j)
	{
	    dh = h[i] - M.h[j];
	    dk = k[i] - M.k[j];
	    nd2[id2++] = dh*dh + dk*dk;
	}
    }
    sort(&nd2[0], &nd2[nd2.size()]);
    double mbadwith = MBadness() + M.MBadness();
    list<int>::iterator idif = ssdIdxFree.begin();
    for (int *pnd2 = &(nd2[0]); pnd2 != &(nd2[nd2.size()]); ++pnd2)
    {
	for ( ; idif != ssdIdxFree.end() && ss->d2hi[*idif] < *pnd2;
		++idif ) {};
	// we found a baddie if we reached the end of ssdIdxFree table
	// or if the current distance is smaller than d2lo
	if (idif == ssdIdxFree.end() || *pnd2 < ss->d2lo[*idif])
	{
	    mbadwith += 2.0;
	}
	// otherwise it is a matching distance
	else
	{
	    ++idif;
	}
    }
    // now add penalty for being outside the SandSphere
    return mbad;
}

////////////////////////////////////////////////////////////////////////
// Molecule operators
////////////////////////////////////////////////////////////////////////

Molecule& Molecule::Shift(double dh, double dk)
{
    subtract_out_penalty();
    for (int i = 0; i < NAtoms; ++i)
    {
	h[i] += (int) round(h[i] + dh);
	k[i] += (int) round(k[i] + dk);
    }
    add_out_penalty();
    set_mbad_abadMax();
    return *this;
}

Molecule& Molecule::Center()
{
    double mean_h = 0.0;
    double mean_k = 0.0;
    for (int i = 0; i < NAtoms; ++i)
    {
	mean_h += h[i];
	mean_k += k[i];
    }
    mean_h /= NAtoms;
    mean_k /= NAtoms;
    Shift(-mean_h, -mean_k);
    return *this;
}

Molecule& Molecule::Rotate(double phi, double h0, double k0)
{
    subtract_out_penalty();
    // define rotation matrix
    double Rz[2][2] = {
	{ cos(phi),	-sin(phi)  },
	{ sin(phi), 	 cos(phi)  }
    };
    double hshift, kshift;
    for (int i = 0; i < NAtoms; ++i)
    {
	// shift coordinates origin to the center of rotation
	double hshift = h[i] - h0;
	double kshift = k[i] - k0;
	// rotate
	double hrot = Rz[0][0]*hshift + Rz[0][1]*kshift;
	double krot = Rz[1][0]*hshift + Rz[1][1]*kshift;
	// shift back
	h[i] = (int) round(hrot + h0);
	k[i] = (int) round(krot + k0);
    }
    add_out_penalty();
    set_mbad_abadMax();
    return *this;
}

Molecule& Molecule::Part(const Molecule& M, const int cidx)
{
    h.resize(1, M.h[cidx]);
    k.resize(1, M.k[cidx]);
    fix_size();
    return *this;
}

Molecule& Molecule::Part(const Molecule& M, const list<int>& cidx)
{
    vector<int> *h_new, *k_new;
    if (&M == this)
    {
	h_new = new vector<int>;
	k_new = new vector<int>;
    }
    else
    {
	h_new = &h;
	k_new = &k;
    }
    h_new->resize(cidx.size());
    k_new->resize(cidx.size());
    int idx = 0;
    for ( list<int>::const_iterator li = cidx.begin();
	    li != cidx.end(); ++li )
    {
	if (*li < 0 || *li >= M.NAtoms)
	{
	    throw range_error("in Molecule::Part()");
	}
	(*h_new)[idx] = M.h[*li];
	(*k_new)[idx] = M.k[*li];
	++idx;
    }
    if (&M == this)
    {
	h = *h_new;
	k = *k_new;
	delete h_new, k_new;
    }
    fix_size();
    return *this;
}

Molecule& Molecule::Pop(const Molecule& M, const int cidx)
{
    if (cidx < 0 || cidx >= M.NAtoms)
    {
	throw range_error("in Molecule::Pop(list<int>)");
    }
    if (this != &M)  *this = M;
    h.erase(h.begin() + cidx);
    k.erase(k.begin() + cidx);
    fix_size();
    return *this;
}

Molecule& Molecule::Pop(const Molecule& M, const list<int>& cidx)
{
    if (cidx.size() == 0)
    {
	if (this != &M)  *this = M;
	return *this;
    }
    list<int> sidx(cidx);
    sidx.sort();
    if (sidx.front() < 0 || sidx.back() >= M.NAtoms)
    {
	throw range_error("in Molecule::Pop(list<int>)");
    }
    sidx.push_back(M.NAtoms);
    vector<int> h_new, k_new;
    int j = 0;
    for ( list<int>::iterator li = sidx.begin();
	    li != sidx.end(); ++li )
    {
	for (; j < *li; ++j)
	{
	    h_new.push_back(M.h[j]);
	    k_new.push_back(M.k[j]);
	}
	j = *li + 1;
    }
    h = h_new;
    k = k_new;
    fix_size();
    return *this;
}

Molecule& Molecule::Clear()
{
    h.clear();
    k.clear();
    fix_size();
    return *this;
}

Molecule& Molecule::Add(Molecule& M)
{
    for (int i = 0; i < M.NAtoms; ++i)
    {
	h.push_back(M.h[i]);
	k.push_back(M.k[i]);
    }
    fix_size();
    return *this;
}

Molecule& Molecule::Add(int nh, int nk)
{
    h.push_back(nh);
    k.push_back(nk);
    fix_size();
    return *this;
}

Molecule& Molecule::MoveAtomTo(int idx, int nh, int nk)
{
    if (idx >= NAtoms)
    {
	throw range_error("in Molecule::MoveAtomTo()");
    }
    UnCache();
    h[idx] = nh;
    k[idx] = nk;
    return *this;
}

struct Molecule::badness_at
{
    badness_at() : h(0), k(0)
    {
	numeric_limits<double> double_info;
	abad = double_info.max();
    }
    badness_at(int nh, int nk, double nbad) :
	h(nh), k(nk), abad(nbad) { }
    int h, k;
    double abad;
};
bool operator<(
	const Molecule::badness_at& lhs, const Molecule::badness_at& rhs
	)
{
    return lhs.abad < rhs.abad;
}

list<Molecule::badness_at> Molecule::find_good_distances(
	int trials, const list<int>& ssd_idx
	)
{
    if (NAtoms > ss->NDist)
    {
	cerr << "E: molecule too large for finding a new position" << endl;
	throw InvalidMolecule();
    }
    // prepare a return list:
    list<badness_at> retlist;
    for (int i = 0; i < trials; ++i)
    {
	int idx = gsl_rng_uniform_int(BGA::rng, ssd_idx.size());
	list<int>::iterator lit = ssd_idx.begin();
	advance(lit, idx);
	double radius = ss->d[*lit];
	double phi = 2.0*M_PI*gsl_rng_uniform(BGA::rng);
	int nh = h[0] + (int)round(radius * cos(phi));
	int nk = k[0] + (int)round(radius * sin(phi));
	retlist.push_back( badness_at(nh, nk, ABadnessAt(nh, nk)) );
    }
    return retlist;
}

list<Molecule::badness_at> Molecule::find_good_triangles(
	int trials, const list<int>& ssd_idx
	)
{
    // try to generate 2 triangles with good distances, this may not always
    // work, so this returns a list, which is either empty, or has even
    // number of entries
    if (NAtoms > ss->NDist)
    {
	cerr << "E: molecule too large for finding a new position" << endl;
	throw InvalidMolecule();
    }
    else if (NAtoms < 2)
    {
	cerr << "E: molecule too small, triangulation not possible" << endl;
	throw InvalidMolecule();
    }
    // prepare a return list:
    list<badness_at> retlist;
    // prepare discrete RNG
    double afit[NAtoms];
    for (int i = 0; i != NAtoms; ++i)
    {
	afit[i] = AFitness(i);
    }
    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(NAtoms, afit);
    for (int nt = 0; nt < trials; ++nt)
    {
	// pick first atom and free distance
	int a1 = gsl_ran_discrete(BGA::rng, table);
	int a2 = gsl_ran_discrete(BGA::rng, table);
	for (int i = 0; i < 5 && a1 == a2; ++i)
	{
	    a2 = gsl_ran_discrete(BGA::rng, table);
	}
	if (a1 == a2)
	{
	    a2 = (a1 + 1) % NAtoms;
	}
	int idf1 = gsl_rng_uniform_int(BGA::rng, ssd_idx.size());
	int idf2 = gsl_rng_uniform_int(BGA::rng, ssd_idx.size()-1) + 1;
	idf2 = (idf2 + idf1) % ssd_idx.size();
	if (idf1 == idf2) throw(runtime_error("idf1 == idf2"));
	list<int>::iterator lit;
	lit = ssd_idx.begin(); advance(lit, idf1);
	double r13 = ss->d[*lit];
	lit = ssd_idx.begin(); advance(lit, idf2);
	double r23 = ss->d[*lit];
	double r12 = dist(a1, a2);
	// is triangle [a1 a2 a3] possible?
	if (r13 + r23 < r12 || fabs(r13 - r23) > r12) continue;
	// here we can construct 2 triangles
	double longdir[2] = {  (h[a2]-h[a1])/r12, (k[a2]-k[a1])/r12 };
	double perpdir[2] = { -longdir[1], 	  longdir[0] };
	double xlong = (r13*r13 + r12*r12 - r23*r23) / (2.0*r12);
	double xperp = sqrt(r13*r13 - xlong*xlong);
	int nh1 = h[a1] + (int) round(xlong*longdir[0] + xperp*perpdir[0]);
	int nk1 = k[a1] + (int) round(xlong*longdir[1] + xperp*perpdir[1]);
	// store the result
	badness_at res1(nh1, nk1, ABadnessAt(nh1, nk1));
	retlist.push_back(res1);
	// jump out if it is a good location or if we did enough trials
	/* debug:
	   cout << nt << " a1 = " << a1 << " a2 = " << a2 << endl;
	   cout << nt << " res1.abad = " << res1.abad << endl;
	   cout << nt << " res1.h = " << res1.h <<
	       " res1.k = " << res1.k << endl;
	   */
	if (res1.abad == 0.0 || ++nt == trials)  break;
	// 2nd triangle:
	int nh2 = h[a1] + (int) round(xlong*longdir[0] - xperp*perpdir[0]);
	int nk2 = k[a1] + (int) round(xlong*longdir[1] - xperp*perpdir[1]);
	// store the result
	badness_at res2(nh2, nk2, ABadnessAt(nh2, nk2));
	retlist.push_back(res2);
	/* debug:
	   cout << nt << " res2.abad = " << res2.abad << endl;
	   cout << nt << " res2.h = " << res2.h <<
	       " res2.k = " << res2.k << endl;
	   cout << nt <<
	   " longdir=[" << longdir[0] << ' ' << longdir[1] << "]," <<
	   " perpdir=[" << perpdir[0] << ' ' << perpdir[1] << "]," <<
	   " xlong=" << xlong << "," <<
	   " xperp=" << xperp << "," <<
	   " r12=" << r12 << "," <<
	   " r13=" << r13 << "," <<
	   " r23=" << r23 << "," <<
	   endl;
	   */
	if (res2.abad == 0.0 || ++nt == trials)  break;
    }
    gsl_ran_discrete_free(table);
    return retlist;
}

Molecule& Molecule::Evolve(int trials)
{
    list<int>::iterator lit;
    if (NAtoms == ss->NAtoms)
    {
	cerr << "E: full-sized molecule cannot Evolve()" << endl;
	throw InvalidMolecule();
    }
    // evolution is trivial for empty or 1-atom molecule
    switch (NAtoms)
    {
	case 0:
	    Add(0, 0);
	    return *this;
	case 1:
	    // make sure ssdIdxFree is updated
	    if (!cached) calc_db();
	    list<badness_at> ba = find_good_distances(1, ssdIdxFree);
	    Add(ba.front().h, ba.front().k);
	    Center();
	    return *this;
    }
    // make sure ssdIdxFree is updated
    if (!cached) calc_db();
    // here we can be sure that NAtoms >= 2
    list<badness_at> trials_log = find_good_distances(trials/3, ssdIdxFree);
    // if we did not find a good place, let's try triangulation
    if (trials_log.back().abad != 0.0)
    {
	list<badness_at> fgt =
	    find_good_triangles(trials-trials_log.size(), ssdIdxFree);
	trials_log.insert(trials_log.end(), fgt.begin(), fgt.end());
    }
    // find the best trial and Add atom to that place
    badness_at *best;
    // maybe the last one is good:
    if (trials_log.back().abad == 0.0)
    {
	best = &trials_log.back();
    }
    else
    {
	best = &(*min_element(trials_log.begin(), trials_log.end()));
    }
    Add(best->h, best->k);
    Center();
    return *this;
}

Molecule& Molecule::Degenerate()
{
    // make sure valarray abad is updated
    if (!cached) calc_db();
    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(NAtoms, &(abad[0]));
    int idx = gsl_ran_discrete(BGA::rng, table);
    gsl_ran_discrete_free(table);
    Pop(idx);
    Center();
    return *this;
}

namespace BGA_Molecule_MateWith
{
    list<int> random_wt_choose(int N, const double *p, int Np)
    {
	list<int> retlst;
	if (N > Np)
	{
	    throw(runtime_error("too many items to choose"));
	}
	// check trivial cases
	else if (N == 0)
	{
	    return retlst;
	}
	else if (N == Np)
	{
	    for (int i = 0; i < Np; ++i)
	    {
		retlst.push_back(i);
	    }
	    return retlst;
	}
	if ( p+Np != find_if(p, p+Np, bind2nd(less<double>(),0.0)) )
	{
	    throw(runtime_error("negative choice probability"));
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
	for (int i = 0, Nprob = Np; i != N; ++i, --Nprob)
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
	    retlst.push_back(val[idx]);
	    // overwrite this element with the last number
	    prob[idx] = prob[Nprob-1];
	    val[idx] = val[Nprob-1];
	}
	return retlst;
    }

    bool compare_MFitness(Molecule& lhs, Molecule& rhs)
    {
	return lhs.MFitness() < rhs.MFitness();
    }
}

Molecule Molecule::MateWith(const Molecule& Male, int trials)
{
    using namespace BGA_Molecule_MateWith;
    if (NAtoms != Male.NAtoms)
    {
	cerr << "E: mating of unequal molecules" << endl;
	throw InvalidMolecule();
    }
    double patom_f[NAtoms];
    for (int i = 0; i < NAtoms; ++i)  patom_f[i] = AFitness(i);
    double patom_m[Male.NAtoms];
    for (int i = 0; i < Male.NAtoms; ++i)  patom_m[i] = Male.AFitness(i);
    list<Molecule> children;
    for (int i = 0; i < trials; ++i)
    {
	Molecule egg(ss);
	list<int> egg_atoms =
	    random_wt_choose(NAtoms, patom_f, NAtoms/2);
	egg.Part(*this, egg_atoms);
	Molecule sperm(ss);
	list<int> sperm_atoms =
	    random_wt_choose(Male.NAtoms, patom_m, NAtoms - egg.NAtoms);
	sperm.Part(Male, sperm_atoms);
	egg.mount(sperm);
	children.push_back(egg);
	if (egg.MBadness() == 0.0)  break;
    }
    // find the best child
    Molecule *best;
    // maybe the last child is a prodigy
    if (children.back().MBadness() == 0.0)
    {
	best = &children.back();
    }
    else
    {
	best = &( *min_element(children.begin(), children.end(),
		compare_MFitness) );
    }
    best->Center();
    return *best;
}

Molecule& Molecule::mount(Molecule& sperm)
{
    if (NAtoms < sperm.NAtoms)
    {
	cerr << "E: sperm molecule is larger than egg" << endl;
	throw InvalidMolecule();
    }
    // mounting is trivial for empty or 1-atom molecules
    if (NAtoms == 0 || sperm.NAtoms ==0)
    {
	return *this;
    }
    else if (NAtoms == 1)
    {
	// here sperm.NAtoms must be equal 1
	// make sure ssdIdxFree is updated
	if (!cached) calc_db();
	list<badness_at> ba = find_good_distances(1, ssdIdxFree);
	Add(ba.front().h, ba.front().k);
	Center();
	return *this;
    }
    // here we can be sure that NAtoms >= 2, sperm.NAtoms >= 1
    // make sure ssdIdxFree is updated for egg and sperm
    if (!cached) calc_db();
    if (!sperm.cached) sperm.calc_db();
    // calculate probabilities of choosing an egg free distance,
    // free distances common to egg and sperm have 5 times higher probability
    const double common_prob = 5.0;
    double peggidx[ssdIdxFree.size()];
    double *pp = peggidx;
    for (list<int>::iterator
	    eidx = ssdIdxFree.begin(), sidx = sperm.ssdIdxFree.begin();
	    eidx != ssdIdxFree.end(); ++eidx, ++pp)
    {
	for (; sidx != sperm.ssdIdxFree.end() && *sidx < *eidx; ++sidx) { }
	// no more common distances if sidx was used
	if (sidx != sperm.ssdIdxFree.end() && *sidx == *eidx)
	{
	    *pp = common_prob;
	}
	else
	{
	    *pp = 1.0;
	}
    }
    // let's start with random mount:
    list<int> d_idx = random_wt_choose(ssdIdxFree.size(), peggidx, d_len);
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
	read_token("NAtoms", NAtoms) &&
	read_token("delta", delta);
    if (!state)
    {
	return;
    }
    if (fmt == "grid")
	format = GRID;
    else if (fmt == "xy")
	format = XY;
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
    // read values to integer vector vhk
    string header;
    vector<int> vhk;
    bool result = read_header(fid, header) && read_data(fid, vhk);
    if (!result) return fid;
    // parse header
    double vhk_scale = 1.0;
    int vhk_NAtoms = vhk.size()/2;
    ParseHeader ph(header);
    if (ph)
    {
	vhk_scale = ph.delta / ss->delta;
	if ( vhk_NAtoms != ph.NAtoms )
	{
	    cerr << "E: " << opened_file << ": expected " << ph.NAtoms <<
		" atoms, read " << vhk_NAtoms << endl;
	    throw IOError();
	}
    }
    // check if all coordinates have been read
    if ( vhk.size() % 2 )
    {
	cerr << "E: " << opened_file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    h.resize(vhk_NAtoms);
    k.resize(vhk_NAtoms);
    for (int i = 0, iv = 0; i < vhk_NAtoms; ++i)
    {
	h[i] = vhk[iv++];
	k[i] = vhk[iv++];
    }
    fix_size();
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

istream& Molecule::ReadXY(istream& fid)
{
    // read values to integer vector vxy
    string header;
    vector<double> vxy;
    bool result = read_header(fid, header) && read_data(fid, vxy);
    if (!result) return fid;
    int vxy_NAtoms = vxy.size()/2;
    // check how many numbers were read
    ParseHeader ph(header);
    if (ph && vxy_NAtoms != ph.NAtoms)
    {
	cerr << "E: " << opened_file << ": expected " << ph.NAtoms <<
	    " atoms, read " << vxy_NAtoms << endl;
	throw IOError();
    }
    if ( vxy.size() % 2 )
    {
	cerr << "E: " << opened_file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    h.resize(vxy.size()/2);
    k.resize(vxy.size()/2);
    for (size_t i = 0, iv = 0; i < vxy.size()/2; ++i)
    {
	h[i] = (int) round(vxy[iv++] / ss->delta);
	k[i] = (int) round(vxy[iv++] / ss->delta);
    }
    fix_size();
    return fid;
}

bool Molecule::ReadXY(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    opened_file = file;
    bool result = ReadXY(fid);
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

bool Molecule::WriteXY(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtXY();
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

Molecule& Molecule::OutFmtXY()
{
    output_format = XY;
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
	case M.XY:
	    result = M.ReadXY(fid);
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
    switch (M.output_format)
    {
	case M.GRID:
	    fid << "# BGA molecule format = grid" << endl;
	    fid << "# NAtoms = " << M.NAtoms << endl;
	    fid << "# delta = " << M.ss->delta << endl;
	    for (int i = 0; i < M.NAtoms; ++i)
	    {
		fid << M.h[i] << '\t' << M.k[i] << endl;
	    }
	    break;
	case M.XY:
	    fid << "# BGA molecule format = xy" << endl;
	    fid << "# NAtoms = " << M.NAtoms << endl;
	    fid << "# delta = " << M.ss->delta << endl;
	    for (int i = 0; i < M.NAtoms; ++i)
	    {
		fid <<
		    M.ss->delta * M.h[i] << '\t' <<
		    M.ss->delta * M.k[i] << endl;
	    }
	    break;
	case M.ATOMEYE:
	    // this format outputs atom badnesses
	    if (!M.cached) M.calc_db();
	    double xyz_lo = 0.0;
	    double xyz_hi = 1.0;
	    double xyz_range = xyz_hi - xyz_lo;
	    if (M.NAtoms > 0)
	    {
		double xyz_extremes[6] = {
		    -M.ss->dmax,
		    1.01*M.ss->delta*(*min_element(M.h.begin(), M.h.end())),
		    1.01*M.ss->delta*(*min_element(M.k.begin(), M.k.end())),
		    M.ss->dmax,
		    1.01*M.ss->delta*(*max_element(M.h.begin(), M.h.end())),
		    1.01*M.ss->delta*(*max_element(M.k.begin(), M.k.end()))
		};
		xyz_lo = *min_element(xyz_extremes, xyz_extremes+6);
		xyz_hi = *max_element(xyz_extremes, xyz_extremes+6);
		xyz_range = xyz_hi - xyz_lo;
	    }
	    fid << "# BGA molecule format = atomeye" << endl;
	    fid << "# NAtoms = " << M.NAtoms << endl;
	    fid << "# delta = " << M.ss->delta << endl;
	    fid << "Number of particles = " << M.NAtoms << endl;
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
	    for (int i = 0; i < M.NAtoms; i++)
	    {
		fid <<
		    (M.h[i]*M.ss->delta - xyz_lo)/xyz_range << " " <<
		    (M.k[i]*M.ss->delta - xyz_lo)/xyz_range << " " <<
		    0.5 << " " <<
		    M.abad[i] << endl;
	    }
	    break;
    }
    return fid;
}

void Molecule::PrintBadness()
{
    // call to MBadness() will update abadMax if necessary
    cout << "MBadness() = " << MBadness() << endl;
    cout << "ABadness() =";
    double mx = (NAtoms > 0) ? abad.max() : 0.0;
    for (int i = 0; i < NAtoms; ++i)
    {
	cout << ' ' << ABadness(i);
	if (ABadness(i) == mx)
	{
	    cout << '*';
	    mx += 1.0;
	}
    }
    cout << endl;
}

void Molecule::PrintFitness()
{
    cout << "MFitness() = " << MFitness() << endl;
    cout << "AFitness() =";
    double mx = (NAtoms > 0) ? abad.max() : 0.0;
    for (int i = 0; i < NAtoms; ++i)
    {
	cout << ' ' << AFitness(i);
	if (ABadness(i) == mx)
	{
	    cout << '*';
	    mx += 1.0;
	}
    }
    cout << endl;
}

////////////////////////////////////////////////////////////////////////
// Population definitions
////////////////////////////////////////////////////////////////////////

void Population::init()
{
    // check whether all members use the same SandSphere
    /*
    iterator ibad =
	adjacent_find(begin(), end(), molecule_SandSpheres_differ);
    if (ibad != end())
    {
	cerr << "E: population contains molecules with different sandsphere" <<
	    endl;
	throw InvalidPopulation();
    }
    */
}

vector<Couple> Population::FindCouples(int NCouples)
{
    double fitness[size()];
    double *pf = fitness;
    for (iterator i = begin(); i != end(); ++i)
    {
	*(pf++) = i->MFitness();
    }
    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(size(), fitness);
    vector<Couple> v(NCouples);
    typedef vector<Couple>::iterator Cit;
    for (Cit i = v.begin(); i != v.end(); ++i)
    {
	i->Male = &at( gsl_ran_discrete(BGA::rng, table) );
	i->Female = &at( gsl_ran_discrete(BGA::rng, table) );
    }
    gsl_ran_discrete_free(table);
    return v;
}
