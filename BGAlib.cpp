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
    d2 = d*d;
    d2lo = d2 - 2*delta*d + delta*delta;
    d2hi = d2 + 2*delta*d + delta*delta;
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

void Molecule::init()
{
    cached = false;
    // check coordinate sizes
    if (h.size() != k.size())
    {
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

typedef struct {
    int d2;
    int i;
    int j;
} d2idx_type;

int d2idx_cmp(const d2idx_type& p, const d2idx_type& q)
{
    if (p.d2 < q.d2)
	return -1;
    else if (p.d2 == q.d2)
	return 0;
    else
	return 1;
}

void Molecule::calc_df()
{
    cached = true;
//    ssdIdxUsed.clear();
    ssdIdxFree.clear();
    d2idx_type d2idx[NDist];
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
    abadMax = abad.max();
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

///* declaration of a pdffit molecule objects */
//class molecule
//{
//public:
//    // constructors
//    molecule(int N = 0); // build empty molecule with Natoms
//    molecule(const char*); // build molecule from file
//    molecule& operator=(const molecule&); // copy assignment
//    // IO functions
//    bool read_stru(const char*); // read molecule from file
//    bool save_stru(const char*); // save molecule to file
//    std::ifstream& operator>>(std::ifstream& s) {return read_stru(s);}
//    friend std::ostream& operator<<(std::ostream& s, molecule& M)
//    {return M.save_stru(s);}
//    std::ostream& save_stru(std::ostream& s) // write molecule to stream
//    {
//	switch(out_fmt) {
//	    case PLAIN:		save_stru_plain(s); break;
//	    case ATOMEYE:	save_stru_aeye(s); break;
//	}
//	return s;
//    }
//    molecule& ofmt_plain() // select plain output format
//    {out_fmt = PLAIN; return *this;} 
//    molecule& ofmt_aeye() // select atom eye output format
//    {out_fmt = ATOMEYE; return *this;}
//    double compare(molecule& m1); // compare 2 molecules:
//    // property interface functions
//    double blen(size_t i, size_t j); // get distance between i-th and j-th atom
//    double blen(size_t k); // get k-th element of the distance matrix
//    double blen_max(); // find maximum distance
//    bool   set_xyz(size_t n, double xn, double yn, double zn); // set position
//    bool   set_all_xyz(std::valarray<double>* pxv, std::valarray<double>* pyv, std::valarray<double>* pzv);
//    inline double get_x(size_t n) { return x[n]; }
//    const  std::valarray<double>& get_x() { return x; };
//    inline double get_y(size_t n) { return y[n]; }
//    const  std::valarray<double>& get_y() { return y; };
//    inline double get_z(size_t n) { return z[n]; }
//    const  std::valarray<double>& get_z() { return z; };
//    void   set_Uiso(size_t n, double Uison) // set temperature factor
//    { Uiso[n] = Uison; }
//    bool   set_all_Uiso(std::valarray<double>* pUisov);
//    void   set_all_Uiso(double Uison) { Uiso = Uison; }
//
//    double get_Uiso(size_t n) {return Uiso[n];}
//    const  std::valarray<double>& get_Uiso() { return Uiso; };
//    inline size_t get_Natoms() {return Natoms;}
//    // other
//    void resize(int Nnew);
//    void   set_verbose(int v) {verbose = v;} // set verbosity level
//    const char* get_title_c_str() {return title.c_str();}
//    mutable std::valarray<size_t> idx_pair; // index of all pairs in Mdist
//    inline size_t ij2idx(size_t i, size_t j) // map matrix to 1D vector
//    {
//	assert(i < Natoms && j < Natoms);
//	return i + j*Natoms;
//    }
//    inline void idx2ij(size_t idx, size_t* i, size_t* j) // map matrix to 1D vector
//    {
//	assert(idx < NMdist);
//	*j = idx / Natoms;
//	*i = idx - *j * Natoms;
//    }
//    std::list<size_t> idx_close_pairs(double rlow); // find joined atoms
//    std::list<size_t> idx_stretched_pairs(double maxblen, double rng); // find lonely atoms
//    std::list<size_t> idx_distant_pairs(double maxdist); // find distant atoms
//    const  std::valarray<double> get_xyzLims() {
//	std::valarray<double> xyzl(xyzLims, 6);
//	return xyzl;
//    }
//private:
//    void   defaults(); // helper constructor function
//    void   alloc(size_t Natoms); // allocate vectors for Natoms
//    double distance(size_t i, size_t j); // calculate distance between i,j
//    // properties
//    double xyzLims[6]; // { xmin, xmax, ymin, ymax, zmin, zmax }
//    std::string title;
//    size_t Natoms; // number of atoms in the unit cell
//    size_t NMdist; // number of atom pairs in the molecule
//    std::vector<std::string> atom_type; // isotope/element of atom i
//    std::valarray<double> x; // coordinate of atom i, 0<=x<1, i<Natoms
//    std::valarray<double> y; // coordinate of atom i, 0<=y<1
//    std::valarray<double> z; // coordinate of atom i, 0<=z<1
//    std::valarray<double> o; // occupancy of atom i
//    std::valarray<double> b; // neutron scattering length for atom i
//    std::valarray<double> Uiso; // isotropic temperature factor in A**2
//    double bavg; // average neutron scattering length
//    bool   orthogonal; // nonzero for orthogonal lattice
//    std::valarray<double> Mdist; // 1D representation of distance matrix
//    std::valarray<bool> good_Mdist; // state of cached Mdist
//    // IO implementation
//    std::string stru_file; // input structure file
//    std::ifstream& read_stru(std::ifstream&); // read molecule in detected format
//    std::ifstream& read_stru_aeye(std::ifstream&); // read structure in atomeye format
//    std::ifstream& read_stru_plain(std::ifstream&); // read structure in plain format
//    std::ifstream& read_stru_pdffit(std::ifstream&); // read structure in PDFFIT format
//    std::ostream& save_stru_plain(std::ostream&); // save in plain format
//    std::ostream& save_stru_aeye(std::ostream&); // save in atomeye CFG format
//    int verbose;
//    // helper functions:
//    template <class T> void grow_shrink(std::valarray<T>& v, size_t n);
//    inline std::string deblank(std::string& s); // removes trailing blanks
//    inline std::string gsub_comma_space(std::string& s); // :s/,/ /g
//    inline void token_read_error(std::ifstream&, std::iostream::pos_type, std::string);
//    void   Error( std::string Message ); // pj: todo error handler
//    enum file_format {PLAIN=1, ATOMEYE};
//    file_format out_fmt;
//};

/***********************************************************************
* Here is what people have been up to:
*
* $Log$
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
