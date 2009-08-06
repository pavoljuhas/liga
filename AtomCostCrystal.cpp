#include "/home/juhas/arch/i686/include/dbprint.h"
/***********************************************************************
 * Short Title: AtomCostCrystal - atom cost calculation with PBC
 *
 * Comments: implementation of class AtomCostCrystal.
 *      AtomCostCrystal inherits from AtomCost.
 *
 * $Id$
 *
 * <license text>
 ***********************************************************************/

#include <cassert>
#include <cmath>
#include "AtomCostCrystal.hpp"
#include "AtomSequence.hpp"
#include "PointsInSphere.hpp"
#include "LigaUtils.hpp"
#include "Lattice.hpp"
#include "Crystal.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// class AtomCostCrystal
////////////////////////////////////////////////////////////////////////

// constructor

AtomCostCrystal::AtomCostCrystal(const Crystal* cluster) : AtomCost(cluster)
{
    use_distances = false;
}

// public methods - overloaded

void AtomCostCrystal::resetFor(const Molecule* clust)
{
    this->AtomCost::resetFor(clust);
    this->arg_cluster = static_cast<const Crystal*>(clust);
    assert(this->arg_cluster == this->AtomCost::arg_cluster);
    assert(!this->use_distances);
    this->_rmax = arg_cluster->getRmax();
    pair<double,double> rext = arg_cluster->getRExtent();
    _sph.reset(new PointsInSphere(rext.first, rext.second,
                arg_cluster->getLattice()) );
}


double AtomCostCrystal::eval(const Atom_t& a, int flags)
{
    return eval(&a, flags);
}


double AtomCostCrystal::eval(const Atom_t* pa, int flags)
{
    // assign arguments
    this->arg_atom = pa;
    this->_gradient_flag = flags & GRADIENT;
    // begin calculation
    resizeArrays();
    resetGradient();
    vector<double>::iterator ptcii = this->partial_costs.begin();
    vector<int>::iterator pcntii = this->pair_counts.begin();
    // reset result data
    this->total_cost = 0.0;
    this->total_pair_count = 0;
    // Cartesian separation vector mapped to unit cell
    R3::Vector rcv;
    for (AtomSequenceIndex seq(arg_cluster); !seq.finished(); seq.next())
    {
	// assertion checks
	assert(ptcii < this->partial_costs.end());
	assert(pcntii < this->pair_counts.end());
        // calculation
        crst_atom = seq.ptr();
        rcv = arg_atom->r - crst_atom->r;
        const pair<double,int> costcount = pairCostCount(rcv);
        *(ptcii++) = costcount.first;
        *(pcntii++) = costcount.second;
        this->total_cost += costcount.first;
        this->total_pair_count += costcount.second;
        bool cutitoff = !this->_gradient_flag && this->apply_cutoff &&
            this->total_cost + arg_atom->Badness() > this->cutoff_cost;
        if (cutitoff)   break;
    }
    bool islowestcost = this->apply_cutoff &&
        arg_atom->Badness() + this->total_cost < this->lowest_cost;
    if (islowestcost)
    {
	this->lowest_cost = arg_atom->Badness() + this->total_cost;
        this->cutoff_cost =
            min(this->cutoff_cost, this->lowest_cost + this->cutoff_range);
    }
    this->_gradient_cached = this->_gradient_flag;
    return this->totalCost();
}


int AtomCostCrystal::totalPairCount() const
{
    return this->total_pair_count;
}


const vector<int>& AtomCostCrystal::pairCounts() const
{
    return this->pair_counts;
}


// public methods - specific

// Evaluate cost of Cartesian separation between two sites
// The separation vector is tranformed to unit cell vector.
// This method also sets

pair<double,int>
AtomCostCrystal::pairCostCount(const R3::Vector& cv, bool skipzero)
{
    const Lattice& lat = arg_cluster->getLattice();
    static R3::Vector ucv;
    ucv = lat.ucvCartesian(cv);
    R3::Vector rc_dd;
    double paircost = 0.0;
    int paircount = 0;
    for (_sph->rewind(); !_sph->finished(); _sph->next())
    {
        rc_dd = ucv + lat.cartesian(_sph->mno());
        double d = R3::norm(rc_dd);
        if (d > this->_rmax)    continue;
        if (skipzero && d < NS_LIGA::eps_distance)  continue;
        double dd = this->nearDistance(d) - d;
        paircost += penalty(dd) * this->getScale();
        paircount += 1;
        if (this->_gradient_flag && d > NS_LIGA::eps_distance)
        {
            static R3::Vector g_dd_xyz;
            g_dd_xyz = (-1.0/d) * rc_dd;
            double g_pcost_dd = penalty_gradient(dd) * this->getScale();
            this->_gradient += g_pcost_dd * g_dd_xyz;
        }
    }
DBPRINT(paircount);
    pair<double,int> rv(paircost, paircount);
    return rv;
}

// protected methods

void AtomCostCrystal::resizeArrays()
{
    size_t sz = arg_cluster->countAtoms();
    bool isresized = (sz == this->partial_costs.size());
    if (isresized)  return;
    this->partial_costs.resize(sz);
    this->pair_counts.resize(sz);
}

// End of file
