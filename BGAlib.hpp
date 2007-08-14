/***********************************************************************
* Short Title: object declarations for single element Liga Algorithm
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
#include <set>
#include "RegisterSVNId.hpp"
#include "Matrix.hpp"
#include "BGAutils.hpp"
#include "R3linalg.hpp"

class AtomCost;

namespace {
RegisterSVNId BGAlib_hpp_id("$Id$");
}

// global random number generator
namespace BGA
{
    extern double eps_badness;
};

// exceptions
struct InvalidDistanceTable { };
struct InvalidMolecule { };
struct InvalidPopulation { };

using namespace std;

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
    iterator find_nearest(const double& d);
    iterator find_nearest_unused(const double& d, valarray<bool>& used);
    iterator return_back(const double&);
    vector<double> unique();
    // data members
    mutable int NAtoms;   	// target number of atoms
    mutable int Nuniqd;   	// number of unique distances
    mutable double max_d;
private:
    void init();
};

enum triangulation_type { LINEAR, PLANAR, SPATIAL, NTGTYPES };

class Atom_t
{
    public:

	Atom_t(double rx0, double ry0, double rz0, double bad0=0.0);
        template <class V> Atom_t(const V& r0, double bad0=0.0);
        R3::Vector r;
	double Badness() const;
	double FreeBadness() const;
	double IncBadness(double db);
	double DecBadness(double db);
	double ResetBadness(double b = 0.0);
	bool fixed;
	triangulation_type ttp;
	int pmxidx;	    // pair matrix index

    private:

	double badness;
};

template <class V>
Atom_t::Atom_t(const V& r0, double bad0) :
    fixed(false), ttp(LINEAR), badness(bad0)
{
    r = r0[0], r0[1], r0[2];
}

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

enum StructureType { MOLECULE, CRYSTAL };

class Molecule
{
    public:

	// friends
	friend class AtomSequence;
	friend class AtomCost;
	friend bool operator==(const Molecule&, const Molecule&);
	friend bool BondAngleFilter_t::Check(Atom_t*, Molecule*);
	friend bool LoneAtomFilter_t::Check(Atom_t*, Molecule*);
	friend ostream& operator<<(ostream& os, Molecule& M);
	friend istream& operator>>(istream& is, Molecule& M);

	// class data
	// fit parameters
	static double tol_dd;
	static double tol_nbad;	// tolerance of normalized badness
	static double tol_r;	// position tolerance in RelaxAtom
	static bool evolve_jump;
	static bool evolve_relax;
	static bool degenerate_relax;
	static double evolve_frac;
	static vector<AtomFilter_t*> atom_filters;
	static double lookout_prob;

	// class methods
	static void OutFmtXYZ();	// output format for operator>>
	static void OutFmtAtomEye();    // output format for operator>>

	// data
	list<int> trace;

	// constructors
	Molecule();
	Molecule(const DistanceTable&);
	Molecule(const DistanceTable&, const int s, const double* px,
		const double* py, const double* pz);
	Molecule(const DistanceTable&, const vector<double>& vx,
		const vector<double>& vy, const vector<double>& vz);
	Molecule(const Molecule& M);		// copy constructor
	Molecule& operator=(const Molecule&);	// assignment

	// destructor
	virtual ~Molecule();

	// methods - class registration and type info
	bool Register();
	virtual StructureType type() const  { return MOLECULE; }
	virtual string typeStr() const { return "molecule"; }

	// methods - fitness/badness evaluation
	inline bool isCloseEnough(const double& dd) const
	{
	    return fabs(dd) < tol_dd;
	}
	double Badness() const;	    // total badness
	double NormBadness() const; // normalized badness
	inline bool Full() const     { return !(NAtoms() < maxNAtoms()); }
	inline int NAtoms() const    { return atoms.size(); }
	inline int maxNAtoms() const { return max_natoms; }
	void setMaxNAtoms(int s);
	inline int NDist()  const
	{
	    int n = NAtoms();
	    return n*(n-1)/2;
	}
	inline double max_dTarget() const { return dTarget.back(); }
	void Recalculate();	    // update everything

	// methods - molecule operations
	void Shift(double dh, double dk, double dl);	// move all atoms
	void Center();	  // center w/r to the center of mass

	// atom operations
	inline const Atom_t& getAtom(const int cidx) { return *(atoms[cidx]); }
	void Pop(const int cidx);	// erase
	void Pop(const list<int>& cidx);
	virtual void Clear();		// remove all atoms
	void Add(const Molecule& M);	// add specified molecule
	void Add(double rx0, double ry0, double rz0);	// add single atom
	void Add(const Atom_t& a);	// add single atom
	void Fix(const int cidx);	// mark atom as fixed
	int NFixed() const;		// count fixed atoms
	void RelaxAtom(const int cidx);	// relax internal atom
	void RelaxAtom(vector<Atom_t*>::iterator);
	void RelaxExternalAtom(Atom_t& a);
	void Evolve(const int* est_triang);
	void Degenerate(int Npop=1);	// Pop Npop atoms with abad[i] weight

	// IO functions
	bool ReadXYZ(const char*); 	// read real coordinates
	bool WriteFile(const char*); 	// save in current output_format
	bool WriteXYZ(const char*); 	// save real coordinates
	bool WriteAtomEye(const char*);	// export in AtomEye format
	void PrintBadness() const;	// total and per-atomic badness
	void PrintFitness();		// total and per-atomic fitness

    protected:

        // data
	DistanceTable dTarget;
	int max_natoms;
	vector<Atom_t*> atoms;		// vector of pointers to atoms

	// methods
	virtual AtomCost* getAtomCostCalculator();
	void addNewAtomPairs(Atom_t* pa);
	void removeAtomPairs(Atom_t* pa);
	int push_good_distances(vector<Atom_t>& vta, double* afit, int ntrials);
	int push_good_triangles(vector<Atom_t>& vta, double* afit, int ntrials);
	int push_good_pyramids(vector<Atom_t>& vta, double* afit, int ntrials);
	int push_second_atoms(vector<Atom_t>& vta, int ntrials);
	int push_third_atoms(vector<Atom_t>& vta, int ntrials);
	valarray<int> good_neighbors_count(const vector<Atom_t>& vta);
	void filter_good_atoms(vector<Atom_t>& vta,
		double evolve_range, double hi_abad);
	bool check_atom_filters(Atom_t*);

    private:

	// types
	enum file_fmt_type {XYZ = 1, ATOMEYE};
	class ParseHeader;

	// class data
	static file_fmt_type output_format;

	// data
	SymmetricMatrix<double> pmx_used_distances;
	SymmetricMatrix<double> pmx_partial_costs;
	set<int> free_pmx_slots;
	mutable double badness;		// molecular badness

	// methods
	// constructor helper
	void init();
	int getPairMatrixIndex();
	// IO helpers
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

#endif	// BGALIB_HPP_INCLUDED
