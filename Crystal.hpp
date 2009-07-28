/***********************************************************************
* Short Title: derived class from Molecule for periodic crystal structure
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef CRYSTAL_T_HPP_INCLUDED
#define CRYSTAL_T_HPP_INCLUDED

#include <boost/shared_ptr.hpp>
#include "Molecule.hpp"

class Lattice;
class AtomCostCrystal;

class Crystal : public Molecule
{
    public:

	// friends
	friend class AtomCostCrystal;

	// data

	// constructors
	Crystal();
	Crystal(const DistanceTable&);
	Crystal(const Crystal&);

        // destructor
	virtual ~Crystal();

        // operators
	virtual Crystal& operator=(const Molecule&);
	Crystal& operator=(const Crystal&);
        virtual Molecule* clone() const;        // create a clone

	// methods - class registration and type info
	virtual StructureType type() const  { return CRYSTAL; }
	virtual std::string typeStr() const { return "crystal"; }

        // methods
        virtual void setDistanceTable(const DistanceTable&);
        virtual void setDistReuse(bool);

        void setLattice(const Lattice&);
        const Lattice& getLattice() const;

        void setRmax(double rmax);
        double getRmax() const;
        std::pair<double,double> getRExtent() const;

        virtual double cost() const;
	virtual int countPairs() const;
        virtual void recalculate() const;
        virtual AtomCost* getAtomCostCalculator() const;

	virtual void Shift(const R3::Vector& drc);

	virtual void Clear();
	virtual void Add(const Atom_t& a);  // add single atom
        virtual const std::pair<int*,int*>& Evolve(const int* est_triang);
	virtual void Degenerate(int Npop=1);

    protected:

        // data
        boost::shared_ptr<DistanceTable> _full_distance_table;
        mutable SymmetricMatrix<int> pmx_pair_counts;

        // methods
	virtual void addNewAtomPairs(Atom_t* pa);
	virtual void removeAtomPairs(Atom_t* pa);
        virtual const TriangulationAnchor&
            getLineAnchor(const RandomWeighedGenerator& rwg);
        virtual const TriangulationAnchor&
            getPlaneAnchor(const RandomWeighedGenerator& rwg);
        virtual const TriangulationAnchor&
            getPyramidAnchor(const RandomWeighedGenerator& rwg);
        virtual void resizePairMatrices(int sz);
        virtual boost::python::object newDiffPyStructure() const;
        virtual void setFromDiffPyStructure(boost::python::object);

    private:

	// crystal specific data
        boost::shared_ptr<const Lattice> _lattice;
	double _rmin;
	double _rmax;
        mutable double _cost;
        mutable int _count_pairs;
        mutable bool _cost_data_cached;
        mutable std::pair<double,double> _rextent;
        mutable bool _rextent_cached;

        // class methods
        boost::shared_ptr<Lattice> getDefaultLattice();

        // methods
        void init();
        void uncacheCostData() const;
        void uncacheRExtent() const;
        void updateRExtent() const;
        void cropDistanceTable();
        const R3::Vector&
            anyOffsetAtomSite(const RandomWeighedGenerator& rwg) const;
        R3::Vector ucvCartesianAdjusted(const R3::Vector& cv) const;
        void shiftToOrigin();
};


#endif	// CRYSTAL_T_HPP_INCLUDED
