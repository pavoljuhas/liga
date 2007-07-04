/***********************************************************************
 * Short Title: AtomCost - functor for external atom cost calculation
 *
 * Comments: implementation of class AtomCost
 *
 * $Id$
 * 
 * <license text>
 ***********************************************************************/

#include <cassert>
#include "AtomCost.hpp"
#include "BGAlib.hpp"
#include "RunPar_t.hpp"

class Molecule;
class RunPar_t;
class Atom_t;

RegisterSVNId AtomCost_cpp_id = "$Id$";

////////////////////////////////////////////////////////////////////////
// class AtomCost
////////////////////////////////////////////////////////////////////////

// constructor

AtomCost::AtomCost(Molecule* m) : arg_atom(NULL)
{
    resetFor(m);
}

// public methods

void AtomCost::resetFor(Molecule* m)
{
    arg_mol = m;
    noCutoff();
    if (arg_mol->maxNAtoms() > int(partial_costs.size()))
    {
	partial_costs.resize(arg_mol->maxNAtoms());
	target_distances.resize(arg_mol->maxNAtoms());
    }
    use_distances = (arg_mol->tol_dd > 0.0);
    if (use_distances && arg_mol->dTarget.size() > useflag.size())
    {
	useflag.resize(arg_mol->dTarget.size(), false);
    }
}

double AtomCost::eval(Atom_t& a)
{
    return eval(&a);
}

double AtomCost::eval(Atom_t* pa)
{
    // assign arguments
    arg_atom = pa;
    // begin calculation
    resetUseFlags();
    total_cost = 0.0;
    vector<double>::iterator tgdii = target_distances.begin();
    vector<double>::iterator ptcii = partial_costs.begin();
    for (AtomSequence seq(arg_mol); !seq.finished(); seq.next())
    {
	// assertion checks
	assert(tgdii < target_distances.end());
	assert(ptcii < partial_costs.end());
	// calculation
	double d = dist(*arg_atom, *seq.ptr());
	vector<double>::iterator pdnear = getNearestDistance(d);
	*(tgdii++) = *pdnear;
	double pcost = penalty(d - *pdnear);
	*(ptcii++) = pcost;
	total_cost += pcost;
	if (apply_cutoff && total_cost + arg_atom->Badness() > cutoff_cost)
	{
	    break;
	}
    }
    if (apply_cutoff && arg_atom->Badness() + total_cost < lowest_cost)
    {
	lowest_cost = arg_atom->Badness() + total_cost;
	cutoff_cost = min(cutoff_cost, lowest_cost + cutoff_range);
    }
    return total_cost;
}

double AtomCost::lowest() const
{
    return lowest_cost;
}

double AtomCost::cutoff() const
{
    return cutoff_cost;
}

void AtomCost::setCutoff(double cf)
{
    apply_cutoff = true;
    cutoff_cost = cf;
    lowest_cost = DOUBLE_MAX;
}

double AtomCost::cutoffRange() const
{
    return cutoff_range;
}

void AtomCost::setCutoffRange(double cutrng)
{
    apply_cutoff = true;
    cutoff_range = cutrng;
    lowest_cost = DOUBLE_MAX;
}

void AtomCost::noCutoff()
{
    apply_cutoff = false;
    lowest_cost = DOUBLE_MAX;
    cutoff_cost = DOUBLE_MAX;
    cutoff_range = DOUBLE_MAX;
}

double AtomCost::total() const
{
    return total_cost;
}

const vector<double>& AtomCost::partialCosts() const
{
    return partial_costs;
}

const vector<double>& AtomCost::targetDistances() const
{
    return target_distances;
}

const vector<int>& AtomCost::usedTargetIndices() const
{
    return useflag_indices;
}

// private methods

void AtomCost::resetUseFlags()
{
    if (useflag_indices.empty())    return;
    for (vector<int>::iterator ii = useflag_indices.begin();
	    ii != useflag_indices.end(); ++ii)
    {
	assert(*ii < int(useflag.size()));
	useflag[*ii] = false;
    }
    useflag_indices.clear();
}

vector<double>::iterator AtomCost::getNearestDistance(const double& d)
{
    DistanceTable& dtgt = arg_mol->dTarget;
    vector<double>::iterator ii = dtgt.find_nearest(d);
    if (!use_distances)
    {
	return ii;
    }
    int idx = ii - dtgt.begin();
    if (useflag[idx])
    {
        int hi, lo, nidx = -1;
        int sz = dtgt.size();
        for (hi = idx + 1; hi < sz && useflag[hi]; ++hi)
        { }
        if (hi < sz)
        {
            nidx = hi;
        }
        for (lo = idx - 1; lo >= 0 && useflag[lo]; --lo)
        { }
        if (lo >= 0 && (nidx < 0 || d - dtgt[lo] < dtgt[nidx] - d))
        {
            nidx = lo;
        }
        idx = nidx;
        ii = dtgt.begin() + nidx;
    }
    if (fabs(*ii - d) < arg_mol->tol_dd)
    {
	useflag[idx] = true;
	useflag_indices.push_back(idx);
    }
    return ii;
}
