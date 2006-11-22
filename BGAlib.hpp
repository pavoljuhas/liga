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
#include <vector>
#include <list>
#include <map>
#include <limits>
#include <gsl/gsl_rng.h>
#include "BGAutils.hpp"

// global random number generator
namespace BGA
{
    extern gsl_rng* rng;
    extern double eps_badness;
};

// exceptions
struct InvalidDistanceTable { };
struct InvalidMolecule { };
struct InvalidPopulation { };

using namespace std;

// constants
const double DOUBLE_MAX = numeric_limits<double>().max();

// helper objects/functions
template <typename T>
struct OrderedPair : pair<T,T>
{
    OrderedPair(const T& x, const T& y) : pair<T,T>(x, y)
    {
	if (pair<T,T>::first > pair<T,T>::second)
	    swap(pair<T,T>::first, pair<T,T>::second);
    }
};
vector<int> random_choose_few(int K, int Np, bool with_repeat = false);
vector<int> random_wt_choose(int K, const double* p, int Np);

inline bool eps_eq(double x, double y)
{
    return fabs(x-y) < BGA::eps_badness;
}

inline bool eps_gt(double x, double y)
{
    return x > y + BGA::eps_badness;
}

inline bool eps_lt(double x, double y)
{
    return x < y - BGA::eps_badness;
}

/* cost function for distance mismatch */
double penalty(double dd);

/* declaration of BGA objects */
class DistanceTable : public vector<double>
{
public:
    // constructors
    DistanceTable();
    DistanceTable(const double*, size_t);
    DistanceTable(const char*);
    DistanceTable(const vector<double>&);
    DistanceTable(const DistanceTable&);
    DistanceTable& operator= (const vector<double>&);
    DistanceTable& operator= (const DistanceTable&);
    // member functions
    vector<double>::iterator find_nearest(const double&);
    vector<double>::iterator return_back(const double&);
    vector<double> unique();
    // data members
    mutable int NAtoms;   	// target number of atoms
    mutable int Nuniqd;   	// number of unique distances
    mutable double max_d;
private:
    void init();
};

class Atom_t
{
public:
    Atom_t(double rx0, double ry0, double rz0, double bad0 = 0.0);
    Atom_t(double r0[3], double bad0 = 0.0);
    Atom_t(valarray<double>& r0, double bad0 = 0.0);
    mutable double r[3];
    double Badness() const;
    double FreeBadness() const;
    double AvgBadness() const;
    double IncBadness(double db);
    double DecBadness(double db);
    double ResetBadness(double b = 0.0);
    bool fixed;
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

class Molecule;

class AtomFilter_t
{
public:
    virtual ~AtomFilter_t() { }
    virtual bool Check(Atom_t*, Molecule* pm = NULL)
    { return true; }
};

class BondAngleFilter_t : public AtomFilter_t
{
public:
    BondAngleFilter_t(double _max_blen) :
	max_blen(_max_blen),
	lo_bangle(0.0),
	hi_bangle(DOUBLE_MAX)
    { }
    virtual ~BondAngleFilter_t() { }
    bool Check(Atom_t*, Molecule* pm);
    double max_blen;
    double lo_bangle;
    double hi_bangle;
};

class LoneAtomFilter_t : public AtomFilter_t
{
public:
    LoneAtomFilter_t(double _max_dist) :
	max_dist(_max_dist)
    { }
    virtual ~LoneAtomFilter_t() { }
    bool Check(Atom_t*, Molecule* pm);
    double max_dist;
};

struct PairDistance_t
{
    void LockTo(Molecule* M, Atom_t* a1, Atom_t* a2);
    void Release(Molecule* M, Atom_t* a1, Atom_t* a2);
    double Badness(Molecule* M, Atom_t* a1, Atom_t* a2);
    double dUsed;
};

class Molecule
{
public:
    // constructors
    Molecule();
    Molecule(const DistanceTable&);
    Molecule(const DistanceTable&, const int s, const double* px,
	    const double* py, const double* pz);
    Molecule(const DistanceTable&, const vector<double>& vx,
	    const vector<double>& vy, const vector<double>& vz);
    Molecule(const Molecule& M);		// copy constructor
    Molecule& operator=(const Molecule&);	// assignment
    ~Molecule();		// destructor
    // fit parameters
    static double tol_dd;
    static double tol_nbad;	// tolerance of normalized badness
    static double tol_r;	// position tolerance in RelaxAtom
    static bool evolve_jump;
    static bool evolve_relax;
    static bool degenerate_relax;
    static double evolve_frac;
    static int center_size;
    static vector<AtomFilter_t*> atom_filters;
    // fitness/badness functions
    double Badness() const;	// total badness
    double NormBadness() const;	// normalized badness
    inline bool Full() const { return !(NAtoms() < max_NAtoms()); }
    // operator functions
    Molecule& Shift(double dh, double dk, double dl);	// move all atoms
    Molecule& Center();	  // center w/r to the center of mass
    // atom operations
    Atom_t Atom(const int cidx);	// get copy of specified atom
    Molecule& Pop(const int cidx);	// erase
    Molecule& Pop(const list<int>& cidx);
    Molecule& Clear();			// remove all atoms
    Molecule& Add(const Molecule& M);	// add specified molecule
    Molecule& Add(double rx0, double ry0, double rz0);	// add single atom
    Molecule& Add(const Atom_t& a);			// add single atom
    Molecule& Fix(const int cidx);		// mark atom as fixed
    int NFixed() const;				// count fixed atoms
    Molecule& RelaxAtom(const int cidx);	// relax internal atom
    Molecule& RelaxAtom(vector<Atom_t*>::iterator);
    void RelaxExternalAtom(Atom_t& a);
    Molecule& Evolve(int ntd1=50, int ntd2=100, int ntd3=500,
	    double lookout_prob=0.0);
    Molecule& Degenerate(int Npop=1);	// Pop Npop atoms with abad[i] weight
    // IO functions
    bool ReadXYZ(const char*); 		// read real coordinates
    bool WriteFile(const char*); 	// save in current output_format
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
    inline int max_NAtoms() const { return val_max_NAtoms; }
    void Set_max_NAtoms(int s);
    inline double max_dTarget() const { return dTarget.back(); }
    // history trace
    list<int> trace;
private:
    // constructor helper
    void init();
    // data storage
    DistanceTable dTarget;
    int val_max_NAtoms;
    // atoms must precede pairs
    vector<Atom_t*> atoms;		      // vector of pointers to atoms
    map<OrderedPair<Atom_t*>,PairDistance_t> pairs;
    friend void PairDistance_t::LockTo(Molecule*, Atom_t*, Atom_t*);
    friend void PairDistance_t::Release(Molecule*, Atom_t*, Atom_t*);
    friend bool operator==(const Molecule&, const Molecule&);
    friend bool BondAngleFilter_t::Check(Atom_t*, Molecule*);
    friend bool LoneAtomFilter_t::Check(Atom_t*, Molecule*);
    // badness evaluation
    mutable double badness;		// molecular badness
    int push_good_distances(vector<Atom_t>& vta, double* afit, int ntrials);
    int push_good_triangles(vector<Atom_t>& vta, double* afit, int ntrials);
    int push_good_pyramids(vector<Atom_t>& vta, double* afit, int ntrials);
    int push_second_atoms(vector<Atom_t>& vta, int ntrials);
    int push_third_atoms(vector<Atom_t>& vta, int ntrials);
    double calc_test_badness(Atom_t& a, double hi_abad = DOUBLE_MAX);
    valarray<int> good_neighbors_count(const vector<Atom_t>& vta);
    void filter_good_atoms(vector<Atom_t>& vta,
	    double evolve_range, double hi_abad);
    bool check_atom_filters(Atom_t*);
    // IO helpers
    enum file_fmt_type {XYZ = 1, ATOMEYE};
    file_fmt_type output_format;
    class ParseHeader;
    istream& ReadXYZ(istream& fid);
    string opened_file;
};
bool operator==(const Molecule& m1, const Molecule& m2);

class Molecule::ParseHeader
{
public:
    ParseHeader(const string& s);
    int NAtoms;
    file_fmt_type format;
    operator bool() {return state;}
private:
    bool state;
    template<typename T> bool read_token(const char* token, T& value);
    const string& header;
};

#endif		// BGALIB_HPP_INCLUDED
