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

        // methods
	Crystal& operator=(const Crystal&);
        void setLattice(const Lattice&);
        const Lattice& getLattice() const;
        void setRRange(double rmin, double rmax);
        std::pair<double,double> getRRange() const;
        std::pair<double,double> getRExtent() const;

    private:

	// crystal specific data
        boost::shared_ptr<Lattice> _lattice;
	double _rmin;
	double _rmax;

        // methods
        void init();

};


#endif	// CRYSTAL_T_HPP_INCLUDED
