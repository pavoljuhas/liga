/***********************************************************************
* Short Title: class Molecule - declaration
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef MOLECULE_HPP_INCLUDED
#define MOLECULE_HPP_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <list>
#include <set>
#include <boost/shared_ptr.hpp>
#include "Atom_t.hpp"
#include "DistanceTable.hpp"
#include "Matrix.hpp"
#include "TraceId_t.hpp"

class AtomFilter_t;
class AtomCost;

enum StructureType { MOLECULE, CRYSTAL };

class Molecule
{
    public:

	// friends
	friend class AtomSequence;
	friend class AtomCost;
	friend bool operator==(const Molecule&, const Molecule&);
	friend class BondAngleFilter_t;
	friend class LoneAtomFilter_t;
	friend std::ostream& operator<<(std::ostream& os, Molecule& M);
	friend std::istream& operator>>(std::istream& is, Molecule& M);

	// class data
	// fit parameters
	static double tol_nbad;	// tolerance of normalized badness
	static double tol_r;	// position tolerance in RelaxAtom
	static bool promotejump;
	static bool promoterelax;
	static bool demoterelax;
	static double promotefrac;
	static std::vector<AtomFilter_t*> atom_filters;
	static double lookout_prob;

	// class methods
	static void OutFmtXYZ();	// output format for operator>>
	static void OutFmtAtomEye();    // output format for operator>>

	// data
        // unique identifier
        const long id;
        std::list<TraceId_t> trace;

	// constructors
	Molecule();
	Molecule(const DistanceTable&);
	Molecule(const DistanceTable&, const int s, const double* px,
		const double* py, const double* pz);
	Molecule(const DistanceTable&,
                const std::vector<double>& vx,
		const std::vector<double>& vy,
                const std::vector<double>& vz);
	Molecule(const Molecule& M);		// copy constructor
	Molecule& operator=(const Molecule&);	// assignment

	// destructor
	virtual ~Molecule();

	// methods - class registration and type info
	bool Register();
	virtual StructureType type() const  { return MOLECULE; }
	virtual std::string typeStr() const { return "molecule"; }

        // methods - molecule configuration
        virtual void setDistanceTable(const DistanceTable&);
        void setDistanceTable(const std::vector<double>&);
        const DistanceTable& getDistanceTable() const;

        virtual void setDistReuse(bool);
        bool getDistReuse() const;

	// methods - fitness/badness evaluation
	virtual double cost() const;    // normalized badness
	const double& Badness() const;	// total badness
        void IncBadness(const double& db) const;
        void DecBadness(const double& db) const;
        void ResetBadness(double b=0.0) const;
	bool full() const;
	int countAtoms() const;
	virtual int countPairs() const;
	int getMaxAtomCount() const;
	void setMaxAtomCount(int cnt);
	double maxTableDistance() const;
        void reassignPairs();	    // improve assignment of distances
	virtual void recalculate() const;   // recalculate everything

	// methods - molecule operations
	void Shift(double dh, double dk, double dl);	// move all atoms
	void Center();	  // center w/r to the center of mass

	// atom operations
	inline const Atom_t& getAtom(const int cidx) { return *(atoms[cidx]); }
	void Pop(const int cidx);	// erase
	void Pop(const std::list<int>& cidx);
	virtual void Clear();		// remove all atoms
	void Add(const Molecule& M);	// add specified molecule
	void Add(double rx0, double ry0, double rz0);	// add single atom
	void Add(const Atom_t& a);	// add single atom
	void Fix(const int cidx);	// mark atom as fixed
	int NFixed() const;		// count fixed atoms
	void RelaxAtom(const int cidx);	// relax internal atom
	void RelaxExternalAtom(Atom_t& a);
        const std::pair<int*,int*>& Evolve(const int* est_triang);
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
        boost::shared_ptr<DistanceTable> _distance_table;
	int _max_atom_count;
        std::vector<Atom_t*> atoms;	// vector of pointers to atoms
	mutable SymmetricMatrix<double> pmx_partial_costs;
	mutable SymmetricMatrix<double> pmx_used_distances;
        mutable std::set<int> free_pmx_slots;
	mutable double _badness;	// molecular badness

	// methods
	virtual AtomCost* getAtomCostCalculator() const;
	virtual void addNewAtomPairs(Atom_t* pa);
	virtual void removeAtomPairs(Atom_t* pa);
        typedef std::vector<Atom_t> AtomArray;
	int push_good_distances(AtomArray& vta, double* afit, int ntrials);
	int push_good_triangles(AtomArray& vta, double* afit, int ntrials);
	int push_good_pyramids(AtomArray& vta, double* afit, int ntrials);
	int push_second_atoms(AtomArray& vta, int ntrials);
	int push_third_atoms(AtomArray& vta, int ntrials);
        std::valarray<int> good_neighbors_count(const AtomArray& vta);
	void filter_good_atoms(AtomArray& vta,
		double evolve_range, double hi_abad);
	bool check_atom_filters(Atom_t*);
        virtual void resizePairMatrices(int sz);

    private:

	// types
	enum file_fmt_type {XYZ = 1, ATOMEYE};
	class ParseHeader;

	// class data
	static file_fmt_type output_format;

        // class methods
        static long getUniqueId();

	// data
        bool _distreuse;

	// methods
	// constructor helper
	void init();
	int getPairMatrixIndex();
        void returnUsedDistances();
	// IO helpers
	std::istream& ReadXYZ(std::istream& fid);
        std::string opened_file;
};

// non-member operators
bool operator==(const Molecule& m1, const Molecule& m2);

class Molecule::ParseHeader
{
public:
    ParseHeader(const std::string& s);
    int NAtoms;
    file_fmt_type format;
    operator bool() {return state;}
private:
    bool state;
    template<typename T> bool read_token(const char* token, T& value);
    const std::string& header;
};

#endif	// MOLECULE_HPP_INCLUDED
