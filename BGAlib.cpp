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

#include "BGAlib.h"

// exceptions
struct InvalidDistanceTable { };
struct IOError { };

////////////////////////////////////////////////////////////////////////
// SandSphere definitions
////////////////////////////////////////////////////////////////////////
SandSphere::SandSphere(int GridMax, const vector<double>& vt) :
    gridmax(GridMax)
{
    init_dist(vt);
}

SandSphere::SandSphere(int GridMax, size_t s, const double *pt) :
    gridmax(GridMax)
{
    vector<double> vt(s);
    for (int i = 0; i < s; ++i) { vt[i] = pt[i]; }
    init_dist(vt);
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
    init_dist(vt);
}


void SandSphere::init_dist(const vector<double>& t)
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
    // fill in and check dist
    dist.resize(NDist);
    for (size_t i = 0; i < NDist; ++i)
    {
	dist[i] = t[i];
    }
    sort(&dist[0], &dist[dist.size()]);
    if (dist[0] <= 0)
    {
	cerr << "E: non-positive entry in target distance table, " <<
	    "dist[0]=" << dist[0] << '\n';
	throw InvalidDistanceTable();
    }
    // calculate grid parameters
    dmax = dist[dist.size()-1];
    delta = dmax/gridmax;
    // calculate d2lo, d2hi
    d2lo.resize(NDist);
    d2hi.resize(NDist);
    d2lo = dist*dist - 2*delta*dist + delta*delta;
    d2hi = dist*dist + 2*delta*dist + delta*delta;
}
 
////////////////////////////////////////////////////////////////////////
// molecule definitions
////////////////////////////////////////////////////////////////////////

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
* Revision 1.1  2005/01/18 23:24:36  juhas
* Initial revision
*
***********************************************************************/
