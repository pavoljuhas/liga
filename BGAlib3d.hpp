/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef BGALIB3D_H_INCLUDED
#define BGALIB3D_H_INCLUDED
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
    const double min_distance = 1.0;
    const int min_distance_penalty = 10;
};

// exceptions
struct InvalidDistanceTable { };
struct InvalidMolecule { };
struct InvalidPopulation { };
struct IOError { };

/* declaration of BGA objects */
using namespace std;

class Molecule;

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
    Atom_t(int h0 = 0, int k0 = 0, int l0 = 0, int bad0 = 0);
    mutable int h, k, l;
    const int Badness();
    const double AvgBadness();
    void IncBadness(int b = 1);
    void DecBadness(int b = 1);
    void ResetBadness();
private:
    int badness;
    int badness_sum;
    int age;
};

struct Pair_t
{
public:
    Pair_t(Molecule *M, Atom_t& a1, Atom_t& a2);
    ~Pair_t();
    // do not allow copying or assignment
    Pair_t(const Pair_t& pair0) {};
    Pair_t& operator=(const Pair_t&) {};
    //
    Atom_t atom1, atom2;
    mutable int d2;
    mutable double d;
private:
    Molecule *pmol;
    int ssdIdxUsed;
    int badness;
};

// Molecule in 2 dimensions
class Molecule
{
public:
    // constructors
    Molecule(SandSphere *SS);
    Molecule(SandSphere *SS, const int s, const int *ph, const int *pk,
	    const int *pl);
    Molecule(SandSphere *SS, const vector<int>& vh, const vector<int>& vk,
	    const vector<int>& vl);
    Molecule(SandSphere *SS, const int s, const double *px, const double *py,
	    const double *pz);
    Molecule(SandSphere *SS,
	    const vector<double>& vx, const vector<double>& vy,
	    const vector<double>& vz);
//    Molecule(const Molecule& M);		// copy constructor
//    Molecule& operator=(const Molecule&);	// assignment
//    ~Molecule();		// destructor
    // parameters
    int NDist;    		// length of distance table
    int NAtoms;   		// current number of atoms
    mutable int MaxAtoms;	// target number of atoms
//    double ABadness(int) const;	// fitness of specified atom
//    double AFitness(int) const;	// badness of specified atom
//    double MBadness() const;	// total badness
//    double MFitness() const;	// total fitness
//    double dist(const int& i, const int& j) const;	// d(i,j)
//    // operator functions
//    Molecule& Shift(double dh, double dk);	// shift all atoms
//    Molecule& Center();			// center w/r to the center of mass
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
//    Molecule& Pop(const Molecule& M, const int cidx);		// get M.Pop()
//    Molecule& Pop(const Molecule& M, const list<int>& cidx);	// get M.Pop()
//    Molecule& Clear();			// remove all atoms
    Molecule& Add(Molecule& M);		// add specified molecule
    Molecule& Add(int nh, int nk, int nl);	// add single atom
    Molecule& Add(Atom_t a);			// add single atom
//    Molecule& MoveAtomTo(int idx, int nh, int nk);	// move 1 atom
//    Molecule& Evolve(int trials = 300);	// Add 1 atom to the right place
//    Molecule& Degenerate(int Npop = 1);	// Pop Npop atoms with abad[i] weight
//    Molecule MateWith(const Molecule& Male, int trials = 500);
//    double ABadnessAt(int nh, int nk) const;  // badness for new atom
//    double MBadnessWith(const Molecule& M) const;   // badness for merging
//    // IO functions
//    bool ReadGrid(const char*); 	// read integer grid coordinates
//    bool WriteGrid(const char*); 	// save integer grid coordinates
//    bool ReadXY(const char*); 		// read real coordinates
//    bool WriteXY(const char*); 		// save real coordinates
//    bool WriteAtomEye(const char*);	// export in AtomEye format
    Molecule& OutFmtGrid();		// output format for operator>>
    Molecule& OutFmtXY();               // output format for operator>>
    Molecule& OutFmtAtomEye();          // output format for operator>>
//    friend ostream& operator<<(ostream& os, Molecule& M);
//    friend istream& operator>>(istream& is, Molecule& M);
//    void PrintBadness();		// total and per-atomic badness
//    void PrintFitness();		// total and per-atomic fitness
    void Recalculate(); 	// update everything
private:
//    // constructor helper
    void init();
//    // data storage
    SandSphere *ss;
    list<Atom_t> atom;			// list of all atoms
    list<Pair_t> pair;			// list of all atom pairs
    mutable list<int> ssdIdxFree;	// available elements in ss.dist
    friend class Pair_t;
//    // badness evaluation
//    mutable double abadMax;		// maximum atom badness
//    mutable double mbadness;		// molecular badness
//public:
//    struct badness_at;
//private:
//    list<badness_at> find_good_distances(int trials, const vector<int>& didx);
//    list<badness_at> find_good_triangles(int trials, const vector<int>& didx);
//    list<badness_at> find_good_triangles2(int trials, const vector<int>& didx);
//    void set_mbad_abadMax() const;
//    // MateWith helpers:
//    Molecule& Molecule::mount(Molecule& Male);
//    // IO helpers
//    enum file_fmt_type {GRID = 1, XY, ATOMEYE};
//    file_fmt_type output_format;
//    class ParseHeader;
//    istream& ReadGrid(istream& fid);
//    istream& ReadXY(istream& fid);
//    string opened_file;
};
//
//class Molecule::ParseHeader
//{
//public:
//    ParseHeader(const string& s);
//    int NAtoms;
//    double delta;
//    file_fmt_type format;
//    operator bool() {return state;}
//private:
//    bool state;
//    template<class T> bool read_token(const char *token, T& value);
//    const string& header;
//};
//
//struct Couple
//{
//    Molecule *Male, *Female;
//};
//
//class Population : public vector<Molecule>
//{
//public:
//    Population() : vector<Molecule>() { init(); }
//    Population(size_type n, SandSphere *SS) :
//	vector<Molecule>(n, Molecule(SS))
//    { init(); }
//    Population(size_type n, const Molecule& M) :
//	vector<Molecule>(n, M)
//    { init(); }
//    Population(Population &P) : vector<Molecule>(P) { init(); }
//    template <class InputIterator>
//	Population(InputIterator first, InputIterator last) :
//	    vector<Molecule>(first, last) { init(); }
//public:
//    vector<Couple> FindCouples(int NCouples = 1);
//private:
//    void init();
//};

#endif
