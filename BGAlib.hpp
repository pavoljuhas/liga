/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef BGALIB_HPP_INCLUDED
#define BGALIB_HPP_INCLUDED

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <valarray>
#include <vector>
#include <list>
#include <map>
#include <gsl/gsl_rng.h>
#include "BGAutils.hpp"

// global random number generator
namespace BGA
{
    extern gsl_rng* rng;
    extern double eps_badness;
    double pow2(double x);
    double well(double x);
};

// exceptions
struct InvalidDistanceTable { };
struct InvalidMolecule { };
struct InvalidPopulation { };

using namespace std;

// helper objects/functions
template <class T>
struct OrderedPair : pair<T,T>
{
    OrderedPair(const T& x, const T& y) : pair<T,T>(x, y)
    { if (!(first<second))  swap(first, second); }
};
list<int> random_choose_few(int K, int Np);
list<int> random_wt_choose(int K, const double *p, int Np);

/* declaration of BGA objects */
class DistanceTable : public vector<double>
{
public:
    // constructors
    DistanceTable();
    DistanceTable(const double*, size_t);
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
{ return sqrt(1.0*dist2(a1, a2)); }

class Molecule;

struct Pair_t
{
public:
    Pair_t(Molecule *M, Atom_t *a1, Atom_t *a2);
    Pair_t(Molecule *M, Atom_t *a1, Atom_t *a2, const Pair_t&);
    ~Pair_t();
    // do not allow copying or assignment
    Pair_t(const Pair_t& pair0);
    Pair_t& operator=(const Pair_t&);
    //
    mutable double d;
private:
    friend class Molecule;
    Atom_t *atom1, *atom2;
    Molecule *owner;
    double dUsed;
    double badness;
};

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
    // fit parameters
    static double tol_dd;
    static double tol_nbad;
    static bool evolve_jump;
    static double evolve_frac;
    static int center_size;
    static double (*penalty)(double);
    // fitness/badness functions
    double Badness() const;	// total badness
    double NormBadness() const;	// normalized badness
    inline bool Full() const { return !(NAtoms() < max_NAtoms()); }
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
    // IO functions
    bool ReadXYZ(const char*); 		// read real coordinates
    bool WriteXYZ(const char*); 	// save real coordinates
    bool WriteAtomEye(const char*);	// export in AtomEye format
    Molecule& OutFmtXYZ();		// output format for operator>>
    Molecule& OutFmtAtomEye();          // output format for operator>>
    friend ostream& operator<<(ostream& os, Molecule& M);
    friend istream& operator>>(istream& is, Molecule& M);
    void PrintBadness();		// total and per-atomic badness
    void PrintFitness();		// total and per-atomic fitness
    void Recalculate(); 	// update everything
    inline int NDist()  const { return pairs.size(); }
    inline int NAtoms() const { return atoms.size(); }
    inline int max_NAtoms() const { return dTarget.NAtoms; }
private:
    // constructor helper
    void init();
    // data storage
    DistanceTable dTarget;
    // atoms must precede pairs
    list<Atom_t> atoms;			// list of all atoms
    map<OrderedPair<Atom_t*>,Pair_t*> pairs;  // map Atom_t* to Pair_t objects
    friend class Pair_t;
    // badness evaluation
    mutable double badness;		// molecular badness
    int push_good_distances(vector<Atom_t>& vta, double *afit, int ntrials);
    int push_good_triangles(vector<Atom_t>& vta, double *afit, int ntrials);
    int push_good_pyramids(vector<Atom_t>& vta, double *afit, int ntrials);
    void calc_test_badness(Atom_t& a);
    void filter_good_atoms(vector<Atom_t>& vta,
	    double evolve_range, double lo_abad);
    // IO helpers
    enum file_fmt_type {XYZ = 1, ATOMEYE};
    file_fmt_type output_format;
    class ParseHeader;
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

#endif		// BGALIB_HPP_INCLUDED
