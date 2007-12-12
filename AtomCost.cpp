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
#include <cmath>
#include "AtomSequence.hpp"
#include "AtomCost.hpp"
#include "Molecule.hpp"
#include "LigaUtils.hpp"
#include "Counter.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////

double penalty(double dd)
{
    static Counter* penalty_calls = Counter::getCounter("penalty_calls");
    penalty_calls->count();
    return dd*dd;
}


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
    arg_cluster = m;
    const DistanceTable& dtgt = arg_cluster->getDistanceTable();
    noCutoff();
    use_distances = !arg_cluster->distreuse;
    if (use_distances && dtgt.size() > useflag.size())
    {
	useflag.resize(dtgt.size(), false);
    }
}

double AtomCost::eval(const Atom_t& a)
{
    return eval(&a);
}

double AtomCost::eval(const Atom_t* pa)
{
    // assign arguments
    arg_atom = pa;
    // begin calculation
    resizeArrays();
    resetUseFlags();
    resetLSQArrays();
    total_cost = 0.0;
    vector<double>::iterator tgdii = target_distances.begin();
    vector<double>::iterator ptcii = partial_costs.begin();
    const DistanceTable& dtgt = arg_cluster->getDistanceTable();
    for (AtomSequenceIndex seq(arg_cluster); !seq.finished(); seq.next())
    {
	// assertion checks
	assert(tgdii < target_distances.end());
	assert(ptcii < partial_costs.end());
	// calculation
	double d = dist(*arg_atom, *seq.ptr());
	size_t nearidx = nearDistanceIndex(d);
	const double& dnear = dtgt[nearidx];
	*(tgdii++) = dnear;
	double dd = dnear - d;
	double pcost = penalty(dd);
	*(ptcii++) = pcost;
	total_cost += pcost;
	if (use_distances)
	{
	    useflag[nearidx] = true;
	    useflag_indices.push_back(nearidx);
	    useatom_indices.push_back(seq.idx());
	}
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

const vector<int>& AtomCost::usedTargetDistanceIndices() const
{
    return useflag_indices;
}

const vector<int>& AtomCost::usedTargetAtomIndices() const
{
    return useatom_indices;
}

size_t AtomCost::lsqComponentsSize() const
{
    if (lsq_anchors.empty())	lsq_anchors = arg_cluster->atoms;
    return lsq_anchors.size();
}

size_t AtomCost::lsqParametersSize() const
{
    return 3;
}

namespace {
inline double sqrt_double(double x) { return sqrt(x); }
}

const vector<double>& AtomCost::lsqWeights() const
{
    if (!lsq_wt.empty())    return lsq_wt;
    lsq_wt.resize(lsqComponentsSize());
    vector<double>::iterator wtii = lsq_wt.begin();
    for (AtomSequence seq(lsq_anchors); !seq.finished(); seq.next())
    {
	// assertion checks
	assert(wtii < lsq_wt.end());
	*(wtii++) = seq.ptr()->Badness();
    }
    lsq_wt = costToFitness(lsq_wt);
    transform(lsq_wt.begin(), lsq_wt.end(), lsq_wt.begin(), sqrt_double);
    return lsq_wt;
}

const vector<double>& AtomCost::lsqComponents() const
{
    if (!lsq_fi.empty())    return lsq_fi;
    lsq_fi.resize(lsqComponentsSize());
    lsq_di.resize(lsqComponentsSize());
    vector<double>::iterator fii = lsq_fi.begin();
    vector<double>::iterator dii = lsq_di.begin();
    vector<double>::const_iterator tgdii = target_distances.begin();
    vector<double>::const_iterator wtii = lsqWeights().begin();
    for (AtomSequence seq(lsq_anchors); !seq.finished();
	    seq.next(), ++fii, ++dii, ++tgdii, ++wtii)
    {
	assert(fii < lsq_fi.end());
	assert(dii < lsq_di.end());
	assert(tgdii < target_distances.end());
	assert(wtii < lsq_wt.end());
	*dii = dist(*arg_atom, *seq.ptr());
	*fii = *wtii * (*tgdii - *dii);
    }
    return lsq_fi;
}

double AtomCost::lsqJacobianGet(size_t m, size_t n) const
{
    // make sure all lsq_vectors are cached
    if (lsq_fi.empty())	    lsqComponents();
    double Jmn;
    Jmn = lsq_wt[m] * (lsq_anchors[m]->r[n] - arg_atom->r[n]) / lsq_di[m];
    return Jmn;
}

// private methods

const vector<Atom_t*>& AtomCost::getClusterAtoms() const
{
    return arg_cluster->atoms;
}

void AtomCost::resizeArrays()
{
    if (arg_cluster->NAtoms() == int(partial_costs.size()))	    return;
    partial_costs.resize(arg_cluster->NAtoms());
    target_distances.resize(arg_cluster->NAtoms());
}

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
    useatom_indices.clear();
}

void AtomCost::resetLSQArrays()
{
    if (lsq_anchors.empty())	return;
    lsq_anchors.clear();
    lsq_di.clear();
    lsq_fi.clear();
    lsq_wt.clear();
}

size_t AtomCost::nearDistanceIndex(const double& d)
{
    const DistanceTable& dtgt = *(arg_cluster->dTarget);
    int idx = dtgt.find_nearest(d) - dtgt.begin();
    if (use_distances && useflag[idx])
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
    }
    return idx;
}

// End of AtomCost.cpp
