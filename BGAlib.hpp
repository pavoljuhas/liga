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

// Molecule in 2 dimensions
class Molecule
{
public:
    // constructors
    Molecule(SandSphere *SS);
    Molecule(SandSphere *SS, int s, int *ph, int *pk);
    Molecule(SandSphere *SS,
	    const vector<int>& vh, const vector<int>& vk);
    Molecule(SandSphere *SS, int s, double *px, double *py);
    Molecule(SandSphere *SS,
	    const vector<double>& vx, const vector<double>& vy);
    Molecule(const Molecule& M);		// copy constructor
    Molecule& operator=(const Molecule&);	// assignment
    ~Molecule();		// destructor
    // parameters
    int NDist;    		// length of distance table
    int NAtoms;   		// current number of atoms
    int MaxAtoms;		// target number of atoms
    double ABadness(int) const;	// fitness of specified atom
    double AFitness(int) const;	// badness of specified atom
    double MBadness() const;	// total badness
    double MFitness() const;	// total fitness
    double dist(const int& i, const int& j) const;	// d(i,j)
    inline int dist2(const int& i, const int& j) const	// squared d(i,j)
    {
	int dh = h[i] - h[j];
	int dk = k[i] - k[j];
	return dh*dh + dk*dk;
    }
    // operator functions
    Molecule& Shift(double dh, double dk);	// shift all atoms
    Molecule& Center();			// center w/r to the center of mass
    Molecule& Rotate(double phi, double h0 = 0.0, double k0 = 0.0);
    template<class T> 
	inline Molecule& Part(T cidx)	// keep only specified atom
	{
	    return Part(*this, cidx);
	}
    // set this Molecule to a Part of M
    Molecule& Part(const Molecule& M, const int cidx); 	
    Molecule& Part(const Molecule& M, const list<int>& cidx);
    template<class T> 
	inline Molecule& Pop(T cidx)	// remove specified atom(s)
	{
	    return Pop(*this, cidx);
	}
    Molecule& Pop(const Molecule& M, const int cidx);		// get M.Pop()
    Molecule& Pop(const Molecule& M, const list<int>& cidx);	// get M.Pop()
    Molecule& Clear();			// remove all atoms
    Molecule& Add(Molecule& M);		// add specified molecule
    Molecule& Add(int nh, int nk);	// add single atom
    Molecule& MoveAtomTo(int idx, int nh, int nk);	// move 1 atom
    Molecule& Evolve(int trials = 300);	// Add 1 atom to the right place
    Molecule& Degenerate();	// Pop 1 atom with abad[i] probability
    Molecule MateWith(const Molecule& Male, int trials = 500);
    double ABadnessAt(int nh, int nk) const;  // badness for new atom
    double MBadnessWith(const Molecule& M) const;   // badness for merging
    // IO functions
    bool ReadGrid(const char*); 	// read integer grid coordinates
    bool WriteGrid(const char*); 	// save integer grid coordinates
    bool ReadXY(const char*); 		// read real coordinates
    bool WriteXY(const char*); 		// save real coordinates
    bool WriteAtomEye(const char*);	// export in AtomEye format
    Molecule& OutFmtGrid();		// output format for operator>>
    Molecule& OutFmtXY();               // output format for operator>>
    Molecule& OutFmtAtomEye();          // output format for operator>>
    friend ostream& operator<<(ostream& os, Molecule& M);
    friend istream& operator>>(istream& is, Molecule& M);
    void PrintBadness();		// total and per-atomic badness
    void PrintFitness();		// total and per-atomic fitness
    // public utility functions
    inline void UnCache() { cached = false; }
private:
    // constructor helper
    void init();
    // data storage
    SandSphere *ss;
    vector<int> h;			// x-coordinates
    vector<int> k;			// y-coordinates
    // badness evaluation
    mutable bool cached;
    mutable valarray<double> abad;	// individual atom badnesses
    mutable double abadMax;		// maximum atom badness
    mutable double mbad;		// molecular badness
    mutable valarray<int> d2;		// sorted table of squared distances 
    double out_penalty(int nh, int nk) const;	// penalty for run-away atoms
//    vector<int> ssdIdxUsed;	// used elements from ss.dist
    mutable list<int> ssdIdxFree;	// available elements in ss.dist
    // operator helper functions
    void fix_size();		// set all sizes consistently with h.size()
    void calc_db() const;	// update distance and fitness tables
public:
    struct badness_at;
private:
    list<badness_at> find_good_distances(int trials, const list<int>& ssdIdx);
    list<badness_at> find_good_triangles(int trials, const list<int>& ssdIdx);
    void subtract_out_penalty() const;
    void add_out_penalty() const;
    void set_mbad_abadMax() const;
    // MateWith helpers:
    Molecule& Molecule::mount(Molecule& Male);
    // IO helpers
    enum file_fmt_type {GRID = 1, XY, ATOMEYE};
    file_fmt_type output_format;
    class ParseHeader;
    istream& ReadGrid(istream& fid);
    istream& ReadXY(istream& fid);
    string opened_file;
};

class Molecule::ParseHeader
{
public:
    ParseHeader(const string& s);
    int NAtoms;
    double delta;
    file_fmt_type format;
    operator bool() {return state;}
private:
    bool state;
    template<class T> bool read_token(const char *token, T& value);
    const string& header;
};

struct Couple
{
    Molecule *Male, *Female;
};

class Population : public vector<Molecule>
{
public:
    Population() : vector<Molecule>() { init(); }
    Population(size_type n, SandSphere *SS) :
	vector<Molecule>(n, Molecule(SS))
    { init(); }
    Population(size_type n, const Molecule& M) :
	vector<Molecule>(n, M)
    { init(); }
    Population(Population &P) : vector<Molecule>(P) { init(); }
    template <class InputIterator>
	Population(InputIterator first, InputIterator last) :
	    vector<Molecule>(first, last) { init(); }
public:
    vector<Couple> FindCouples(int NCouples = 1);
private:
    void init();
};

#endif
