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

void AtomCostCrystal::resetFor(const Crystal* clust)
{
    arg_cluster = clust;
    noCutoff();
    pair<double,double> r_range = arg_cluster->getRRange();
    this->_rmin = r_range.first;
    this->_rmax = r_range.second;
    pair<double,double> rext = arg_cluster->getRExtent();
    _sph.reset(new PointsInSphere(rext.first, rext.second,
                arg_cluster->getLattice()) );
    this->_lsq_component_size = -1;
}

double AtomCostCrystal::eval(const Atom_t* pa)
{
    // assign arguments
    arg_atom = pa;
    // begin calculation
    resizeArrays();
    total_cost = 0.0;
    const Lattice& lat = arg_cluster->getLattice();
    // find corresponding atom site in the unit cell
    R3::Vector aucv = lat.ucvCartesian(arg_atom->r);
    _lsq_component_size = 0;
    for (AtomSequenceIndex seq(arg_cluster); !seq.finished(); seq.next())
    {
        const Atom_t& a = seq.ref();
        R3::Vector rc_dd;
        for (_sph->rewind(); !_sph->finished(); _sph->next())
        {
            rc_dd = a.r + lat.cartesian(_sph->mno) - aucv;
            double d = R3::norm(rc_dd);
            if (d < this->_rmin || d > this->_rmax || d == 0.0)     continue;
            _lsq_component_size++;
            double dd = this->nearDistance(d) - d;
            double pcost = penalty(dd);
            assert(seq.idx() < int(partial_costs.size()));
            partial_costs[seq.idx()] += pcost;
            total_cost += pcost;
        }
        bool cutitoff = this->apply_cutoff &&
            this->total_cost + arg_atom->Badness() > this->cutoff_cost;
        if (cutitoff)
        {
            _lsq_component_size = -1;
            break;
        }
    }
    bool isnewlowestcost = this->apply_cutoff &&
        arg_atom->Badness() + this->total_cost < this->lowest_cost;
    if (isnewlowestcost)
    {
	this->lowest_cost = arg_atom->Badness() + this->total_cost;
        this->cutoff_cost =
            min(this->cutoff_cost, this->lowest_cost + this->cutoff_range);
    }
    return this->total_cost;
}

size_t AtomCostCrystal::lsqComponentsSize() const
{
    if (_lsq_component_size < 0)
    {
        const char* emsg = "Undetermined number of LSQ components.";
        throw runtime_error(emsg);
    }
    return size_t(_lsq_component_size);
}

const vector<double>& AtomCostCrystal::lsqComponents() const
{
    if (!lsq_fi.empty())    return lsq_fi;
    lsq_fi.resize(lsqComponentsSize());
    lsq_di.resize(lsqComponentsSize());
    lsq_wt.resize(lsqComponentsSize());
    _lsq_jacobian.resize(lsqComponentsSize(), lsqParametersSize());
    const Lattice& lat = arg_cluster->getLattice();
    // find corresponding atom site in the unit cell
    R3::Vector aucv = lat.ucvCartesian(arg_atom->r);
    int lsqidx = 0;
    vector<double>::iterator fii = lsq_fi.begin();
    vector<double>::iterator dii = lsq_di.begin();
    vector<double>::iterator wtii = lsq_wt.begin();
    for (AtomSequenceIndex seq(arg_cluster); !seq.finished(); seq.next())
    {
        const Atom_t& a = seq.ref();
        R3::Vector rc_dd;
        double wt = sqrt(costToFitness(a.Badness()));
        for (_sph->rewind(); !_sph->finished(); _sph->next())
        {
            rc_dd = a.r + lat.cartesian(_sph->mno) - aucv;
            double d = R3::norm(rc_dd);
            if (d < this->_rmin || d > this->_rmax || d == 0.0)   continue;
            lsqidx++;
            assert(fii < lsq_fi.end());
            assert(dii < lsq_di.end());
            assert(wtii < lsq_wt.end());
            *(dii++) = d;
            *(wtii++) = wt;
            double dd = this->nearDistance(d) - d;
            *(fii++) = wt * dd;
            for (int j = 0; j < int(lsqParametersSize()); ++j)
            {
                _lsq_jacobian(lsqidx, j) = wt * (rc_dd[j]) / d;
            }
        }
    }
    return lsq_fi;
}

double AtomCostCrystal::lsqJacobianGet(size_t m, size_t n) const
{
    // make sure all lsq_vectors are cached
    if (lsq_fi.empty())	    lsqComponents();
    return _lsq_jacobian(int(m), int(n));
}

// public methods - specific

double AtomCostCrystal::costOfLattice() const
{
    double totcost = 0.0;
    size_t npts = 0;
    for (_sph->rewind(); !_sph->finished(); _sph->next())
    {
        double d = _sph->r();
        if (d < this->_rmin || d > this->_rmax || d == 0.0)   continue;
        double dd = this->nearDistance(d) - d;
        totcost += penalty(dd);
        npts++;
    }
    double latcost = npts ? totcost/npts : 0.0;
    return latcost;
}


// protected methods

void AtomCostCrystal::resizeArrays()
{
    partial_costs.resize(arg_cluster->NAtoms());
    fill(partial_costs.begin(), partial_costs.end(), 0.0);
}

// End of file
