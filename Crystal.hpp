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
        virtual const DistanceTable& getDistanceTable() const;

        void setLattice(const Lattice&);
        const Lattice& getLattice() const;

        void setRRange(double rmin, double rmax);
        std::pair<double,double> getRRange() const;
        std::pair<double,double> getRExtent() const;

        double cost() const;
        double costOfLattice() const;

    protected:

        // methods
        virtual AtomCost* getAtomCostCalculator() const;

    private:

	// crystal specific data
        boost::shared_ptr<Lattice> _lattice;
        boost::shared_ptr<DistanceTable> _distance_table;
	double _rmin;
	double _rmax;
        mutable double _cost;
        mutable bool _cost_cached;
        mutable double _cost_of_lattice;
        mutable double _cost_of_lattice_cached;

        // methods
        void init();
        void uncacheCost();

};


#endif	// CRYSTAL_T_HPP_INCLUDED
