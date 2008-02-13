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
        virtual void resetFor(const Molecule* clust);
	virtual double eval(const Atom_t* pa);
	int totalPairCount() const;
	const std::vector<int>& pairCounts() const;
	virtual size_t lsqComponentsSize() const;
	virtual const std::vector<double>& lsqComponents() const;
	virtual double lsqJacobianGet(size_t m, size_t n) const;

        // public methods - specific
        std::pair<double,int>
            pairCostCount(const R3::Vector& cv, bool skipzero=false) const;

    protected:

	// data - results
        int total_pair_count;
        std::vector<int> pair_counts;

        /* delete once working
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
        R3::Vector arg_rcuc;    // cartesian positions offset to unit cell

        // data - for intermediate cost evaluation
        // maximum r for PDF range
        double _rmax;
        // lattice points sequencer
        std::auto_ptr<PointsInSphere> _sph;
        int _lsq_component_size;
        mutable blitz::Array<double,2> _lsq_jacobian;
        mutable int _count_evaluated_pairs;

	// private methods
	void resizeArrays();

};  // class AtomCostCrystal

#endif	// ATOMCOSTCRYSTAL_HPP_INCLUDED
