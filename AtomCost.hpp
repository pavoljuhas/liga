/***********************************************************************
* Short Title: AtomCost - functor for external atom cost calculation
*
* Comments: base class for AtomCostCrystal
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef ATOMCOST_HPP_INCLUDED
#define ATOMCOST_HPP_INCLUDED

#include <vector>
#include "RegisterSVNId.hpp"

class Molecule;
class Atom_t;

namespace {
RegisterSVNId AtomCost_hpp_id = "$Id$";
}

class AtomCost
{
    public:

	// data - arguments
	Molecule* arg_mol;
	const Atom_t* arg_atom;

	// constructor
	AtomCost(Molecule* m);

	// destructor
	virtual ~AtomCost() { }

	// public methods
	virtual void resetFor(Molecule* m);
	double eval(const Atom_t& pa);
	virtual double eval(const Atom_t* pa);
	double lowest() const;
	double cutoff() const;
	void setCutoff(double cf);
	double cutoffRange() const;
	void setCutoffRange(double cutrng);
	void noCutoff();
	double total() const;
	const std::vector<double>& partialCosts() const;
	const std::vector<double>& targetDistances() const;
	const std::vector<int>& usedTargetDistanceIndices() const;
	const std::vector<int>& usedTargetAtomIndices() const;
	virtual size_t lsqComponentsSize() const;
	size_t lsqParametersSize() const;
	const std::vector<double>& lsqWeights() const;
	const std::vector<double>& lsqComponents() const;
	double lsqJacobianGet(size_t m, size_t n) const;

    protected:

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

	// protected methods

    private:

	// data members
	std::vector<bool> useflag;

	// private methods
	void resizeArrays();
	void resetUseFlags();
	void resetLSQArrays();
	size_t nearDistanceIndex(const double& d);

};  // class AtomCost

#endif	// ATOMCOST_HPP_INCLUDED
