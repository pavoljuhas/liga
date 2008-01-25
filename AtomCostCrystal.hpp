/***********************************************************************
* Short Title: AtomCostCrystal - functor for cost calculation with PBC.
*
* Comments: 
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef ATOMCOSTCRYSTAL_HPP_INCLUDED
#define ATOMCOSTCRYSTAL_HPP_INCLUDED

#include <vector>
#include <blitz/array.h>
#include "AtomCost.hpp"
#include "R3linalg.hpp"
#include "PointsInSphere.hpp"

class Crystal;
class Atom_t;

class AtomCostCrystal : public AtomCost
{
    public:

	// constructor
	AtomCostCrystal(const Crystal* m);

	// destructor
	virtual ~AtomCostCrystal() { }

	// public methods - overloaded
        virtual void resetFor(const Crystal* clust);
	virtual double eval(const Atom_t* pa);
	virtual size_t lsqComponentsSize() const;
	virtual const std::vector<double>& lsqComponents() const;
	virtual double lsqJacobianGet(size_t m, size_t n) const;

        // public methods - specific
        double costOfLattice() const;

    protected:

        /* delete once working
	// data - results
	bool use_distances;
	bool apply_cutoff;
	double lowest_cost;
	double cutoff_cost;
	double cutoff_range;
	double total_cost;
	std::vector<double> partial_costs;
	std::vector<double> target_distances;
	std::vector<int> useflag_indices;
	std::vector<int> useatom_indices;
	// LSQ specific data
	mutable std::vector<Atom_t*> lsq_anchors;   // anchor atoms
	mutable std::vector<double> lsq_di;	    // model distances
	mutable std::vector<double> lsq_fi;	    // LSQ components
	mutable std::vector<double> lsq_wt;	    // weights
        */

	// protected methods

    private:

	// data - arguments
	const Crystal* arg_cluster;
	const Atom_t* arg_atom;
        R3::Vector arg_rcuc;    // cartesian positions offset to unit cell

        // data - for intermediate cost evaluation
        // PDF range, actual and extended for sphere summation
        double _rmin;
        double _rmax;
        std::auto_ptr<PointsInSphere> _sph;
        int _lsq_component_size;
        mutable blitz::Array<double,2> _lsq_jacobian;

	// private methods
	void resizeArrays();

};  // class AtomCostCrystal

#endif	// ATOMCOSTCRYSTAL_HPP_INCLUDED
