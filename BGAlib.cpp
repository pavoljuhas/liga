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

#include "BGAlib.hpp"

// exceptions
struct InvalidDistanceTable { };
struct InvalidMolecule { };
struct IOError { };

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
	cerr << "E: unable to open `" << file << "'\n";
	throw IOError();
    }
    vector<double> vt;
    // ignore header lines up to the first number:
    string line;
    istringstream istrs;
    double x;
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
		vt.push_back(x);
	    } while (istrs >> x);
	}
    }
    // read the rest of the file
    while (fid >> x)
    {
	vt.push_back(x);
    }
    // check if everything was read
    if (!fid.eof())
    {
	cerr << "E: " << file << ':' << fid.tellg() << ": invalid number\n";
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
	cerr << "E: target distance table is empty\n";
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
    setGridTol(defGridTol);
}

void SandSphere::setGridTol(double t)
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
    UnCache();
    fixsize();
}

void Molecule::fixsize()
{
    UnCache();
    if (h.size() != k.size())
    {
	cerr << "E: invalide coordinate vectors\n";
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
	cerr << "E: molecule too large\n";
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
    // now add penalty for outside SandSphere
    double Ri;
    for (int i = 0; i < NAtoms; ++i)
    {
	Ri = sqrt(h[i]*h[i] + k[i]*k[i] + 0.0);
	if (Ri > ss->gridmax)
	{
	    abad[i] += Ri - ss->gridmax;
	}
    }
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
    static double mbadCache;
    if (!cached)
    {
	calc_df();
	mbadCache = abad.sum();
    }
    return mbadCache;
}

double Molecule::MFitness()
{
    // this will update abadMax if necessary
    double mbadness = MBadness();
    return NAtoms*abadMax - mbadness;
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

////////////////////////////////////////////////////////////////////////
// Molecule operators
////////////////////////////////////////////////////////////////////////
Molecule& Molecule::Shift(int dh, int dk)
{
    for (int i = 0; i < NAtoms; ++i)
    {
	h[i] += dh;
	k[i] += dk;
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
    fixsize();
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
    fixsize();
    return *this;
}

Molecule& Molecule::Pop(const int cidx)
{
    h.erase(h.begin() + cidx);
    k.erase(k.begin() + cidx);
    fixsize();
    return *this;
}


/***********************************************************************
* Here is what people have been up to:
*
* $Log$
* Revision 1.11  2005/01/26 05:41:01  juhas
* added Molecule operations: Shift(int dx, int dy), Center(),
*     Part(list<int>), Pop(list<int>), Pop(int)
* Molecule resize operations grouped to fixsize()
*
* Revision 1.10  2005/01/25 23:19:11  juhas
* added functions for dealing with GridTol,
* SandSphere keeps the list of molecules using the same grid
*
* Revision 1.9  2005/01/25 20:12:29  juhas
* fixed evaluation of d2lo, d2hi, abadMax
*
* Revision 1.8  2005/01/25 17:23:24  juhas
* added Molecule constructors from real coordinates
*
* Revision 1.7  2005/01/25 16:42:33  juhas
* *** empty log message ***
*
* Revision 1.6  2005/01/25 16:23:07  juhas
* simplified comparison function d2idx_cmp
*
* Revision 1.5  2005/01/25 16:16:19  juhas
* BGAlib.h renamed to BGAlib.hpp
*
* Revision 1.4  2005/01/25 04:46:39  juhas
* added penalty for atoms outside SandSphere
*
* Revision 1.3  2005/01/25 04:38:15  juhas
* afit renamed to abad, afitMax renamed to abadMax
* added Molecule definitions for
*     calc_df(),  ABadness(),  AFitness(),  MBadness(),
*     MFitness(),  dist(),  dist2()
*
* Revision 1.2  2005/01/19 21:19:38  juhas
* init_dist() renamed to init()
* list "dist" renamed to "d"
* added a couple of Molecule definitions
*
* Revision 1.1.1.1  2005/01/18 23:24:36  juhas
* BGA - Biosphere Genetic Algorithm
*
***********************************************************************/
