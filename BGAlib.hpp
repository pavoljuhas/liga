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
#include <gsl/gsl_rng.h>

// global random number generator
namespace BGA
{
    extern gsl_rng* rng;
};

// exceptions
struct InvalidDistanceTable { };
struct InvalidMolecule { };
struct InvalidPopulation { };
struct IOError { };

/* declaration of BGA objects */
using namespace std;

class Molecule;

class DistanceTable : public vector<double>
{
public:
    // constructors
    DistanceTable();
    DistanceTable(const double*, size_t);
    template <class InputIterator>
	DistanceTable(InputIterator first, InputIterator last);
    DistanceTable(const char*);
    DistanceTable(const vector<double>&);
    DistanceTable& operator= (const vector<double>&);
    // member functions
    vector<double>::iterator find_nearest(const double&);
    vector<double>::iterator return_back(const double&);
    // data members
    mutable int NAtoms;   	// target number of atoms
    mutable double max_d;
private:
    void init();
};

class SandSphere
{
public:
    // constructors
    SandSphere(int GridMax, const vector<double>& t);
    SandSphere(int GridMax, int s, const double *t);
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
    void SetGridTol(double t);	// set grid tolerance
    double GridTol();
private:
    static const double defGridTol = 1.0;
    double vGridTol;
    // helper functions
    void init(const vector<double>& t);
    list<Molecule *> molecules;
    friend class Molecule;
};

class Atom_t
{
public:
    Atom_t(double rx0, double ry0, double rz0, double bad0 = 0);
    mutable double rx, ry, rz;
    double Badness() const;
    double AvgBadness() const;
    double IncBadness(double db = 1.0);
    double DecBadness(double db = 1.0);
    double ResetBadness(double b = 0.0);
private:
    double badness;
    double badness_sum;
    int age;
};

bool operator==(const Atom_t& a1, const Atom_t& a2);
double dist2(const Atom_t& a1, const Atom_t& a2);
inline double dist(const Atom_t& a1, const Atom_t& a2)
{
    return sqrt(1.0*dist2(a1, a2));
}

struct Pair_t
{
public:
    Pair_t(Molecule *M, Atom_t& a1, Atom_t& a2);
    ~Pair_t();
    // do not allow copying or assignment
    Pair_t(const Pair_t& pair0);
    Pair_t& operator=(const Pair_t&);
    //
    mutable double d2;
    mutable double d;
    Atom_t *atom1, *atom2;
private:
    Molecule *owner;
    double dUsed;
    double badness;
};

// Molecule in 2 dimensions
class Molecule
{
public:
    // constructors
    Molecule(const DistanceTable&);
    Molecule(const DistanceTable&, const int s, const double *px,
	    const double *py, const double *pz);
    Molecule(const DistanceTable&, const vector<double>& vx,
	    const vector<double>& vy, const vector<double>& vz);
    Molecule(const Molecule& M);		// copy constructor
    Molecule& operator=(const Molecule&);	// assignment
    ~Molecule();		// destructor
    // parameters
//    int ABadness(int) const;	// fitness of specified atom
//    int AFitness(int) const;	// badness of specified atom
    double Badness() const;	// total badness
    double Fitness() const;	// total fitness
    double MaxABadness() const;	// total fitness
    // operator functions
    Molecule& Shift(double dh, double dk, double dl);	// move all atoms
    Molecule& Center();	  // center w/r to the center of mass
//    Molecule& Rotate(double phi, double h0 = 0.0, double k0 = 0.0);
//    template<class T> 
//	inline Molecule& Part(T cidx)	// keep only specified atom
//	{
//	    return Part(*this, cidx);
//	}
//    // set this Molecule to a Part of M
//    Molecule& Part(const Molecule& M, const int cidx); 	
//    Molecule& Part(const Molecule& M, const list<int>& cidx);
//    template<class T> 
//	inline Molecule& Pop(T cidx)	// remove specified atom(s)
//	{
//	    return Pop(*this, cidx);
//	}
    // remove atom(s)
    Molecule& Pop(list<Atom_t>::iterator ai);
    Molecule& Pop(const int cidx);
    Molecule& Pop(const list<int>& cidx);
    Molecule& Clear();			// remove all atoms
    Molecule& Add(Molecule& M);		// add specified molecule
    Molecule& Add(double rx0, double ry0, double rz0);	// add single atom
    Molecule& Add(Atom_t a);				// add single atom
//    Molecule& MoveAtomTo(int idx, int nh, int nk);	// move 1 atom
    Molecule& Evolve(int ntd1=50, int ntd2=100, int ntd3=5);
    Molecule& Degenerate(int Npop=1);	// Pop Npop atoms with abad[i] weight
//    Molecule MateWith(const Molecule& Male, int trials = 500);
//    double ABadnessAt(int nh, int nk) const;  // badness for new atom
//    double MBadnessWith(const Molecule& M) const;   // badness for merging
    // IO functions
    bool ReadGrid(const char*); 	// read integer grid coordinates
    bool WriteGrid(const char*); 	// save integer grid coordinates
    bool ReadXYZ(const char*); 		// read real coordinates
    bool WriteXYZ(const char*); 	// save real coordinates
    bool WriteAtomEye(const char*);	// export in AtomEye format
    Molecule& OutFmtGrid();		// output format for operator>>
    Molecule& OutFmtXYZ();		// output format for operator>>
    Molecule& OutFmtAtomEye();          // output format for operator>>
    friend ostream& operator<<(ostream& os, Molecule& M);
    friend istream& operator>>(istream& is, Molecule& M);
    void PrintBadness();		// total and per-atomic badness
    void PrintFitness();		// total and per-atomic fitness
    void Recalculate(); 	// update everything
    inline int NDist()  const { return pairs.size(); }
    inline int NAtoms() const { return atoms.size(); }
private:
    // constructor helper
    void init();
    // data storage
    DistanceTable dTarget;
    //pj:remove ss
    SandSphere *ss;
    // atoms must precede pairs
    list<Atom_t> atoms;			// list of all atoms
    list<Pair_t*> pairs;		// list of all atom Pair_t objects
    //pj: remove ssdIdxFree
    mutable vector<int> ssdIdxFree;	// available elements in ss.dist
    friend class Pair_t;
    mutable int max_NAtoms;		// target number of atoms
    // badness evaluation
    mutable double max_abad;		// maximum atom badness
    mutable double badness;		// molecular badness
    int push_good_distances(vector<Atom_t>& vta, double *afit, int ntrials);
    int push_good_triangles(vector<Atom_t>& vta, double *afit, int ntrials);
    int push_good_pyramids(vector<Atom_t>& vta, double *afit, int ntrials);
    void calc_test_badness(Atom_t& a);
//    list<badness_at> find_good_distances(int trials, const vector<int>& didx);
//    list<badness_at> find_good_triangles(int trials, const vector<int>& didx);
//    list<badness_at> find_good_triangles2(int trials, const vector<int>& didx);
//    void set_mbad_abadMax() const;
//    // MateWith helpers:
//    Molecule& Molecule::mount(Molecule& Male);
    // IO helpers
    enum file_fmt_type {GRID = 1, XYZ, ATOMEYE};
    file_fmt_type output_format;
    class ParseHeader;
    istream& ReadGrid(istream& fid);
    istream& ReadXYZ(istream& fid);
    string opened_file;
};

class Molecule::ParseHeader
{
public:
    ParseHeader(const string& s);
    int NAtoms;
    file_fmt_type format;
    operator bool() {return state;}
private:
    bool state;
    template<class T> bool read_token(const char *token, T& value);
    const string& header;
};


template<class T> typename list<T>::iterator list_at(const list<T>& lst, int n);
list<int> random_choose_few(int K, int Np);
list<int> random_wt_choose(int K, const double *p, int Np);
double vdnorm(const valarray<double>&);
double vddot(const valarray<double>&, const valarray<double>&);
valarray<double> vdcross(const valarray<double>&, const valarray<double>&);

#endif
