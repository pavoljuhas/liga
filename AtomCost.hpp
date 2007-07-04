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
class RunPar_t;
class Atom_t;

namespace {
RegisterSVNId AtomCost_hpp_id = "$Id$";
}

class AtomCost
{
    public:

	// data members

	// constructor
	AtomCost(Molecule* m);

	// destructor
	virtual ~AtomCost() { }

	// public methods
	virtual void resetFor(Molecule* m);
	double eval(Atom_t& pa);
	virtual double eval(Atom_t* pa);
	double lowest() const;
	double cutoff() const;
	void setCutoff(double cf);
	double cutoffRange() const;
	void setCutoffRange(double cutrng);
	void noCutoff();
	double total() const;
	const std::vector<double>& partialCosts() const;
	const std::vector<double>& targetDistances() const;
	const std::vector<int>& usedTargetIndices() const;

    protected:

	// data - arguments
	Molecule* arg_mol;
	Atom_t* arg_atom;
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

	// protected methods

    private:

	// data members
	std::vector<bool> useflag;

	// private methods
	void resetUseFlags();
	std::vector<double>::iterator getNearestDistance(const double& d);

};  // class AtomCost

#endif	// ATOMCOST_HPP_INCLUDED
