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
        const double& getRmax() const;
        std::pair<double,double> getRExtent() const;

        virtual double cost() const;
	virtual int countPairs() const;
        virtual void recalculate() const;

	virtual void Clear();
	virtual void Add(const Atom_t& a);  // add single atom

    protected:

        // data
        mutable SymmetricMatrix<int> pmx_pair_counts;

        // methods
        virtual AtomCost* getAtomCostCalculator() const;
	virtual void addNewAtomPairs(Atom_t* pa);
	virtual void removeAtomPairs(Atom_t* pa);
        virtual const TriangulationAnchor&
            getLineAnchor(const RandomWeighedGenerator& rwg);
        virtual const TriangulationAnchor&
            getPlaneAnchor(const RandomWeighedGenerator& rwg);
        virtual const TriangulationAnchor&
            getPyramidAnchor(const RandomWeighedGenerator& rwg);
        virtual void resizePairMatrices(int sz);

    private:

	// crystal specific data
        boost::shared_ptr<Lattice> _lattice;
	double _rmin;
	double _rmax;
        mutable double _cost;
        mutable int _count_pairs;
        mutable bool _cost_data_cached;

        // methods
        void init();
        void uncacheCostData();
        const R3::Vector&
            anyOffsetAtomSite(const RandomWeighedGenerator& rwg) const;

};


#endif	// CRYSTAL_T_HPP_INCLUDED
