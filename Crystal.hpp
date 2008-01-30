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
	virtual ~Crystal();

        // operators
	Crystal& operator=(const Crystal&);

        // methods
        virtual void setDistanceTable(const DistanceTable&);
        virtual void setDistReuse(bool);

        void setLattice(const Lattice&);
        const Lattice& getLattice() const;

        void setRRange(double rmin, double rmax);
        std::pair<double,double> getRRange() const;
        std::pair<double,double> getRExtent() const;

        virtual double cost() const;
	virtual int countPairs() const;
        virtual void recalculate() const;

	virtual void Clear();

    protected:

        // data
        mutable SymmetricMatrix<int> pmx_pair_counts;

        // methods
        virtual AtomCost* getAtomCostCalculator() const;
	virtual void addNewAtomPairs(Atom_t* pa);
	virtual void removeAtomPairs(Atom_t* pa);
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

};


#endif	// CRYSTAL_T_HPP_INCLUDED
