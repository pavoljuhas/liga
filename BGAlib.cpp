/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/


//#include <iostream>
//#include <string>
//#include <sstream>
//#include <fstream>
//#include <valarray>
//#include <vector>
//#include <list>

#include <stdexcept>
#include "BGAlib.hpp"

// exceptions
struct InvalidDistanceTable { };
struct InvalidMolecule { };
struct IOError { };

// read vector from stream after skipping header lines
template<class T> istream& read_header(istream& fid, vector<T>& v)
{
    // prepare v
    v.clear();
    // ignore header lines up to the first number:
    string line;
    istringstream istrs;
    T x;
    while (!fid.eof() && getline(fid, line))
    {
	istrs.clear();
	istrs.str(line);
	if (! (istrs >> x) )
	{
	    continue;
	}
	else
	{
	    do {
		v.push_back(x);
	    } while (istrs >> x);
	}
    }
    // read the rest of the file
    while (fid >> x)
    {
	v.push_back(x);
    }
    return fid;
}

////////////////////////////////////////////////////////////////////////
// SandSphere definitions
////////////////////////////////////////////////////////////////////////
SandSphere::SandSphere(int GridMax, const vector<double>& vt) :
    gridmax(GridMax)
{
    init(vt);
}

SandSphere::SandSphere(int GridMax, size_t s, const double *pt) :
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
    if (!fid) {
	cerr << "E: unable to read `" << file << "'" << endl;
	throw IOError();
    }
    // read values to vt
    vector<double> vt;
    read_header(fid, vt);
    // check if everything was read
    if (!fid.eof())
    {
	cerr << "E: " << file << ':' << fid.tellg() <<
	    ": invalid number" << endl;
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
    for (size_t i = 0; i < NDist; ++i)
    {
	d[i] = t[i];
    }
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
    h.clear();
    k.clear();
    init();
}

Molecule::Molecule(SandSphere *SS,
	size_t s, int *ph, int *pk) : ss(SS)
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
	size_t s, double *px, double *py) : ss(SS)
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
    for (int i = 0; i < h.size(); ++i)
    {
	h[i] = (int) round(vx[i] / ss->delta);
	k[i] = (int) round(vy[i] / ss->delta);
    }
    init();
}

void Molecule::init()
{
    ss->molecules.push_back(this);
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
    abad.resize(NAtoms, 0.0);
    abadMax = 0.0;
    d2.resize(NDist, 0);
//    ssdIdxUsed.clear();
    ssdIdxFree.clear();
}

Molecule::~Molecule()
{
    ss->molecules.remove(this);
    // debug: cout << "ss->molecules.size() = " << ss->molecules.size() << endl;
}

double Molecule::dist(const int& i, const int& j)
{
    return sqrt(1.0*dist2(i, j));
}

inline int Molecule::dist2(const int& i, const int& j)
{
    int dh = h[i] - h[j];
    int dk = k[i] - k[j];
    return dh*dh + dk*dk;
}

double Molecule::out_penalty(int i)
{
    double Rout = sqrt(h[i]*h[i] + k[i]*k[i] + 0.0) - ss->gridmax;
    return (Rout > 0.0) ? Rout : 0.0;
}

typedef struct {
    int d2;
    int i;
    int j;
} d2idx_type;

bool d2idx_cmp(const d2idx_type& p, const d2idx_type& q)
{
    return (p.d2 < q.d2);
}

void Molecule::calc_df()
{
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
    sort(d2idx, d2idx+NDist, d2idx_cmp);
    for (ij = 0; ij < NDist; ++ij)
    {
	d2[ij] = d2idx[ij].d2;
    }
    // evaluate abad[i]
    abad = 0.0;
    int ssdIdx = 0;
    for (ij = 0; ij < NDist; ++ij)
    {
	for(; ss->d2hi[ssdIdx] < d2idx[ij].d2 && ssdIdx < ss->NDist; ++ssdIdx)
	{
	    ssdIdxFree.push_back(ssdIdx);
	}
	// did we find a baddie?
	if (ssdIdx < ss->NDist  &&  ss->d2lo[ssdIdx] > d2idx[ij].d2)
	{
	    abad[d2idx[ij].i]++;
	    abad[d2idx[ij].j]++;
	}
	// otherwise it is a matching distance
	ssdIdx++;
    }
    // now add penalty for being outside the SandSphere
    for (int i = 0; i < NAtoms; ++i)
    {
	abad[i] += out_penalty(i);
    }
    mbad = abad.sum();
    abadMax = max(abad.max(), (double) NAtoms);
}

double Molecule::ABadness(int i)
{
    if (!cached)
    {
	calc_df();
    }
    return abad[i];
}

double Molecule::AFitness(int i)
{
    // this will update abadMax if necessary
    double badness_i = ABadness(i);
    return abadMax - badness_i;
}

double Molecule::MBadness()
{
    if (!cached)
    {
	calc_df();
    }
    return mbad;
}

double Molecule::MFitness()
{
    // this will update abadMax if necessary
    double mbadness = MBadness();
    return NAtoms*abadMax - mbadness;
}

////////////////////////////////////////////////////////////////////////
// Molecule operators
////////////////////////////////////////////////////////////////////////
Molecule& Molecule::Shift(int dh, int dk)
{
    for (int i = 0; i < NAtoms; ++i)
    {
	if (cached)
	{
	    abad[i] -= out_penalty(i);
	}
	h[i] += dh;
	k[i] += dk;
	if (cached)
	{
	    abad[i] += out_penalty(i);
	}
    }
    if (cached)
    {
	abadMax = max(abad.max(), (double) NAtoms);
	mbad = abad.sum();
    }
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
    Shift( (int) round(-mean_h), (int) round(-mean_k) );
    return *this;
}

Molecule& Molecule::Part(const list<int>& cidx)
{
    vector<int> h_new, k_new;
    for ( list<int>::const_iterator li = cidx.begin();
	    li != cidx.end(); ++li )
    {
	h_new.push_back(h[*li]);
	k_new.push_back(k[*li]);
    }
    h = h_new;
    k = k_new;
    fix_size();
    return *this;
}

Molecule& Molecule::Pop(const list<int>& cidx)
{
    list<int> sidx(cidx);
    sidx.sort();
    sidx.push_back(NAtoms);
    vector<int> h_new, k_new;
    int j = 0;
    for ( list<int>::iterator li = sidx.begin();
	    li != sidx.end(); ++li )
    {
	for (; j < *li; ++j)
	{
	    h_new.push_back(h[j]);
	    k_new.push_back(k[j]);
	}
	j = *li + 1;
    }
    h = h_new;
    k = k_new;
    fix_size();
    return *this;
}

Molecule& Molecule::Pop(const int cidx)
{
    h.erase(h.begin() + cidx);
    k.erase(k.begin() + cidx);
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

Molecule& Molecule::Add(Molecule& m)
{
    for (int i = 0; i < m.NAtoms; ++i)
    {
	h.push_back(m.h[i]);
	k.push_back(m.k[i]);
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

Molecule& Molecule::MoveAtom(int idx, int nh, int nk)
{
    if (idx > NAtoms)
    {
	throw range_error("in Molecule::MoveAtom()");
    }
    UnCache();
    h[idx] = nh;
    k[idx] = nk;
    return *this;
}

////////////////////////////////////////////////////////////////////////
// Molecule IO functions
////////////////////////////////////////////////////////////////////////

bool Molecule::ReadGrid(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid) {
	cerr << "E: unable to read `" << file << "'" << endl;
	throw IOError();
    }
    // read values to integer vector vhk
    vector<int> vhk;
    bool result = read_header(fid, vhk);
    fid.close();
    // check how many numbers were read
    if ( vhk.size() % 2 )
    {
	cerr << "E: " << file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    h.resize(vhk.size()/2);
    k.resize(vhk.size()/2);
    for (int i = 0, iv = 0; i < vhk.size()/2; ++i)
    {
	h[i] = vhk[iv++];
	k[i] = vhk[iv++];
    }
    fix_size();
    return result;
}

bool Molecule::ReadXY(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid) {
	cerr << "E: unable to read `" << file << "'" << endl;
	throw IOError();
    }
    // read values to integer vector vxy
    vector<double> vxy;
    bool result = read_header(fid, vxy);
    fid.close();
    // check how many numbers were read
    if ( vxy.size() % 2 )
    {
	cerr << "E: " << file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    h.resize(vxy.size()/2);
    k.resize(vxy.size()/2);
    for (int i = 0, iv = 0; i < vxy.size()/2; ++i)
    {
	h[i] = (int) round(vxy[iv++] / ss->delta);
	k[i] = (int) round(vxy[iv++] / ss->delta);
    }
    fix_size();
    return result;
}

bool write_file(const char* file, Molecule& m)
{
    // open file for writing
    ofstream fid(file);
    if (!fid) {
	cerr << "E: unable to write to `" << file << "'" << endl;
	throw IOError();
    }
    bool result = (fid << m);
    fid.close();
    return result;
}

bool Molecule::WriteGrid(const char* file)
{
    out_fmt_type org_ofmt = output_format;
    OutFmtGrid();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
    return result;
}

bool Molecule::WriteXY(const char* file)
{
    out_fmt_type org_ofmt = output_format;
    OutFmtXY();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
    return result;
}

bool Molecule::WriteAtomEye(const char* file)
{
    out_fmt_type org_ofmt = output_format;
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

ostream& operator<<(ostream& os, Molecule& m)
{
    switch (m.output_format)
    {
	case m.GRID:
	    os << "# BGA molecule format: grid" << endl;
	    os << "# NAtoms = " << m.NAtoms << endl;
	    os << "# delta = " << m.ss->delta << endl;
	    for (int i = 0; i < m.NAtoms; ++i)
	    {
		os << m.h[i] << '\t' << m.k[i] << endl;
	    }
	    break;
	case m.XY:
	    os << "# BGA molecule format: xy" << endl;
	    os << "# NAtoms = " << m.NAtoms << endl;
	    os << "# delta = " << m.ss->delta << endl;
	    for (int i = 0; i < m.NAtoms; ++i)
	    {
		os <<
		    m.ss->delta * m.h[i] << '\t' <<
		    m.ss->delta * m.k[i] << endl;
	    }
	    break;
	case m.ATOMEYE:
	    os << "# BGA molecule format: aeye" << endl;
	    os << "# NAtoms = " << m.NAtoms << endl;
	    os << "# delta = " << m.ss->delta << endl;
	    os << "Number of particles = " << m.NAtoms << endl;
	    os << "A = 1.0 Angstrom (basic length-scale)" << endl;
	    os << "H0(1,1) = " << 2.0 * m.ss->dmax << " A" << endl;
	    os << "H0(1,2) = 0 A" << endl;
	    os << "H0(1,3) = 0 A" << endl;
	    os << "H0(2,1) = 0 A" << endl;
	    os << "H0(2,2) = " << 2.0 * m.ss->dmax << " A" << endl;
	    os << "H0(2,3) = 0 A" << endl;
	    os << "H0(3,1) = 0 A" << endl;
	    os << "H0(3,2) = 0 A" << endl;
	    os << "H0(3,3) = " << 2.0 * m.ss->dmax << " A" << endl;
	    os << ".NO_VELOCITY." << endl;
	    // 4 entries: x, y, z, Uiso
	    os << "entry_count = 4" << endl;
	    os << "auxiliary[0] = abad [au]" << endl;
	    os << endl;
	    // pj: now it only works for a single Carbon atom in the molecule
	    os << "12.0111" << endl;
	    os << "C" << endl;
	    for (int i = 0; i < m.NAtoms; i++)
	    {
		os <<
		    (m.h[i]*m.ss->delta + m.ss->dmax)/(2*m.ss->dmax) << " " <<
		    (m.k[i]*m.ss->delta + m.ss->dmax)/(2*m.ss->dmax) << " " <<
		    0.5 << " " <<
		    m.abad[i] << endl;
	    }
	    break;
    }
    return os;
}
