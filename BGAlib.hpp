/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef BGALIB_H_INCLUDED
#define BGALIB_H_INCLUDED
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <valarray>
#include <vector>
#include <list>

/* declaration of BGA objects */
using namespace std;

class Molecule;

class SandSphere
{
public:
    // constructors
    SandSphere(int GridMax, const vector<double>& t);
    SandSphere(int GridMax, size_t s, const double *t);
    SandSphere(int GridMax, const char *file);
    // grid parameters
    double delta;		// grid step
    double dmax;		// maximum target distance
    int gridmax;		// maximum grid coordinate
    // target distance table
    int NDist;    		// length of distance table
    int NAtoms;   		// target number of atoms
    valarray<double> d;		// sorted list of target distances
    valarray<double> d2;	// sorted list of squared target distances
    valarray<double> d2lo;	// low limits for d2
    valarray<double> d2hi;	// high limits for d2
    void setGridTol(double t);	// set grid tolerance
    double GridTol();
private:
    static const double defGridTol = 1.0;
    double vGridTol;
    // helper functions
    void init(const vector<double>& t);
    list<Molecule *> molecules;
    friend class Molecule;
};

// Molecule in 2 dimensions
class Molecule
{
public:
    // constructors
    Molecule(SandSphere *SS);
    Molecule(SandSphere *SS, size_t s, int *ph, int *pk);
    Molecule(SandSphere *SS,
	    const vector<int>& vh, const vector<int>& vk);
    Molecule(SandSphere *SS, size_t s, double *px, double *py);
    Molecule(SandSphere *SS,
	    const vector<double>& vx, const vector<double>& vy);
    // parameters
    int NDist;    		// length of distance table
    int NAtoms;   		// target number of atoms
    double ABadness(int);	// fitness of specified atom
    double AFitness(int);	// badness of specified atom
    double MBadness();		// total badness
    double MFitness();		// total fitness
    double dist(const int& i, const int& j);		// d(i,j)
    inline int dist2(const int& i, const int& j);	// squared d(i,j)
    // operator functions
    Molecule& Shift(int dh, int dk);	// shift all atoms
    Molecule& Center();			// center w/r to the center of mass
    Molecule& Part(const list<int>& cidx);	// keep only specified atoms
    Molecule& Pop(const list<int>& cidx);	// remove specified atoms
    Molecule& Pop(const int cidx);	// remove 1 atom
    Molecule& Add(Molecule& m);		// add specified molecule
    Molecule& Add(int nh, int nk);	// add single atom
    Molecule& MoveAtom(int idx, int nh, int nk);	// move 1 atom
    // IO functions
    bool Read(const char*); 	// read integer grid coordinates
    bool Save(const char*); 	// save integer grid coordinates
    bool SaveAeye(const char*); // export coordinates in AtomEye format
    // public utility functions
    inline void UnCache() { cached = false; }
    ~Molecule();		// destructor
private:
    // data storage
    SandSphere *ss;
    vector<int> h;		// x-coordinates
    vector<int> k;		// y-coordinates
    bool cached;
    valarray<double> abad;	// individual atom badnesses
    double abadMax;		// maximum atom badness
    valarray<int> d2;		// sorted table of squared distances 
//    vector<int> ssdIdxUsed;	// used elements from ss.dist
    list<int> ssdIdxFree;	// available elements in ss.dist
    // helper functions
    void init();		// constructor helper
    void calc_df();		// update distance and fitness tables
};

/*
    std::ifstream& operator>>(std::ifstream& s) {return read_stru(s);}
    friend std::ostream& operator<<(std::ostream& s, molecule& M)
    {return M.save_stru(s);}
    std::ostream& save_stru(std::ostream& s) // write molecule to stream
    {
	switch(out_fmt) {
	    case PLAIN:		save_stru_plain(s); break;
	    case ATOMEYE:	save_stru_aeye(s); break;
	}
	return s;
    }
    molecule& ofmt_plain() // select plain output format
    {out_fmt = PLAIN; return *this;} 
    molecule& ofmt_aeye() // select atom eye output format
    {out_fmt = ATOMEYE; return *this;}
    double compare(molecule& m1); // compare 2 molecules:
    // property interface functions
    double blen(size_t i, size_t j); // get distance between i-th and j-th atom
    double blen(size_t k); // get k-th element of the distance matrix
    double blen_max(); // find maximum distance
    bool   set_xyz(size_t n, double xn, double yn, double zn); // set position
    bool   set_all_xyz(std::valarray<double>* pxv, std::valarray<double>* pyv, std::valarray<double>* pzv);
    inline double get_x(size_t n) { return x[n]; }
    const  std::valarray<double>& get_x() { return x; };
    inline double get_y(size_t n) { return y[n]; }
    const  std::valarray<double>& get_y() { return y; };
    inline double get_z(size_t n) { return z[n]; }
    const  std::valarray<double>& get_z() { return z; };
    void   set_Uiso(size_t n, double Uison) // set temperature factor
    { Uiso[n] = Uison; }
    bool   set_all_Uiso(std::valarray<double>* pUisov);
    void   set_all_Uiso(double Uison) { Uiso = Uison; }

    double get_Uiso(size_t n) {return Uiso[n];}
    const  std::valarray<double>& get_Uiso() { return Uiso; };
    inline size_t get_Natoms() {return Natoms;}
    // other
    void resize(int Nnew);
    void   set_verbose(int v) {verbose = v;} // set verbosity level
    const char* get_title_c_str() {return title.c_str();}
    mutable std::valarray<size_t> idx_pair; // index of all pairs in Mdist
    inline size_t ij2idx(size_t i, size_t j) // map matrix to 1D vector
    {
	assert(i < Natoms && j < Natoms);
	return i + j*Natoms;
    }
    inline void idx2ij(size_t idx, size_t* i, size_t* j) // map matrix to 1D vector
    {
	assert(idx < NMdist);
	*j = idx / Natoms;
	*i = idx - *j * Natoms;
    }
    std::list<size_t> idx_close_pairs(double rlow); // find joined atoms
    std::list<size_t> idx_stretched_pairs(double maxblen, double rng); // find lonely atoms
    std::list<size_t> idx_distant_pairs(double maxdist); // find distant atoms
    const  std::valarray<double> get_xyzLims() {
	std::valarray<double> xyzl(xyzLims, 6);
	return xyzl;
    }
private:
    void   defaults(); // helper constructor function
    void   alloc(size_t Natoms); // allocate vectors for Natoms
    double distance(size_t i, size_t j); // calculate distance between i,j
    // properties
    double xyzLims[6]; // { xmin, xmax, ymin, ymax, zmin, zmax }
    std::string title;
    size_t Natoms; // number of atoms in the unit cell
    size_t NMdist; // number of atom pairs in the molecule
    std::vector<std::string> atom_type; // isotope/element of atom i
    std::valarray<double> x; // coordinate of atom i, 0<=x<1, i<Natoms
    std::valarray<double> y; // coordinate of atom i, 0<=y<1
    std::valarray<double> z; // coordinate of atom i, 0<=z<1
    std::valarray<double> o; // occupancy of atom i
    std::valarray<double> b; // neutron scattering length for atom i
    std::valarray<double> Uiso; // isotropic temperature factor in A**2
    double bavg; // average neutron scattering length
    bool   orthogonal; // nonzero for orthogonal lattice
    std::valarray<double> Mdist; // 1D representation of distance matrix
    std::valarray<bool> good_Mdist; // state of cached Mdist
    // IO implementation
    std::string stru_file; // input structure file
    std::ifstream& read_stru(std::ifstream&); // read molecule in detected format
    std::ifstream& read_stru_aeye(std::ifstream&); // read structure in atomeye format
    std::ifstream& read_stru_plain(std::ifstream&); // read structure in plain format
    std::ifstream& read_stru_pdffit(std::ifstream&); // read structure in PDFFIT format
    std::ostream& save_stru_plain(std::ostream&); // save in plain format
    std::ostream& save_stru_aeye(std::ostream&); // save in atomeye CFG format
    int verbose;
    // helper functions:
    template <class T> void grow_shrink(std::valarray<T>& v, size_t n);
    inline std::string deblank(std::string& s); // removes trailing blanks
    inline std::string gsub_comma_space(std::string& s); // :s/,/ /g
    inline void token_read_error(std::ifstream&, std::iostream::pos_type, std::string);
    void   Error( std::string Message ); // pj: todo error handler
    enum file_format {PLAIN=1, ATOMEYE};
    file_format out_fmt;
};
*/

#endif

/***********************************************************************
* Here is what people have been up to:
*
* $Log$
* Revision 1.8  2005/01/25 23:19:11  juhas
* added functions for dealing with GridTol,
* SandSphere keeps the list of molecules using the same grid
*
* Revision 1.6  2005/01/25 17:23:24  juhas
* added Molecule constructors from real coordinates
*
* Revision 1.5  2005/01/25 16:42:33  juhas
* *** empty log message ***
*
* Revision 1.4  2005/01/25 04:39:50  juhas
* afit renamed to abad, afitMax renamed to abadMax
* declared Molecule::calc_df()
*
* Revision 1.3  2005/01/19 21:20:39  juhas
* init_dist() renamed to init()
* list "dist" renamed to "d"
* added a couple of Molecule definitions
*
* Revision 1.2  2005/01/19 00:16:15  juhas
* added some Molecule declarations
*
* Revision 1.1.1.1  2005/01/18 23:24:36  juhas
* BGA - Biosphere Genetic Algorithm
*
*
***********************************************************************/
