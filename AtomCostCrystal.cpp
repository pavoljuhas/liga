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
#include "AtomSequence.hpp"
#include "AtomCostCrystal.hpp"
#include "PointsInSphere.hpp"
#include "LigaUtils.hpp"
#include "Lattice.hpp"
#include "Crystal.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// class AtomCostCrystal
////////////////////////////////////////////////////////////////////////

// constructor

AtomCostCrystal::AtomCostCrystal(Crystal* cluster) : AtomCost(cluster)
{
    use_distances = false;
}

// public methods - overloaded

void AtomCostCrystal::resetFor(Crystal* clust)
{
    arg_cluster = clust;
    noCutoff();
    cluster_atoms = getClusterAtoms();
    cluster_atoms.push_back(NULL);
    pair<double,double> r_range = arg_cluster->getRRange();
    rmin = r_range.first;
    rmax = r_range.second;
    pair<double,double> r_extent = arg_cluster->getRExtent();
    sph.reset(new PointsInSphere(r_extent.first, r_extent.second,
                arg_cluster->getLattice()) );
    _lsq_component_size = -1;
}

double AtomCostCrystal::eval(const Atom_t* pa)
{
    // assign arguments
    arg_atom = pa;
    // begin calculation
    resizeArrays();
    total_cost = 0.0;
    const DistanceTable& dtgt = arg_cluster->getDistanceTable();
    const Lattice& lat = arg_cluster->getLattice();
    // add copy of argument atom to calculate self contribution
    Atom_t arg_atom_uc = shiftToUnitCell(*arg_atom);
    cluster_atoms.back() = &arg_atom_uc;
    _lsq_component_size = 0;
    for (AtomSequenceIndex seq(cluster_atoms); !seq.finished(); seq.next())
    {
        R3::Vector rc_dd;
        for (sph->rewind(); !sph->finished(); sph->next())
        {
            rc_dd = seq.ptr()->r + lat.cartesian(sph->mno) - arg_atom_uc.r;
            double d = R3::norm(rc_dd);
            if (d < rmin || d > rmax || d == 0.0)   continue;
            _lsq_component_size++;
            size_t nearidx = nearDistanceIndex(d);
            const double& dnear = dtgt[nearidx];
            double dd = dnear - d;
            double pcost = penalty(dd);
            assert(seq.idx() < int(partial_costs.size()));
            partial_costs[seq.idx()] += pcost;
            total_cost += pcost;
        }
        if (apply_cutoff && total_cost + arg_atom->Badness() > cutoff_cost)
        {
            _lsq_component_size = -1;
            break;
        }
    }
    if (apply_cutoff && arg_atom->Badness() + total_cost < lowest_cost)
    {
	lowest_cost = arg_atom->Badness() + total_cost;
	cutoff_cost = min(cutoff_cost, lowest_cost + cutoff_range);
    }
    cluster_atoms.back() = NULL;
    return total_cost;
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
    const DistanceTable& dtgt = arg_cluster->getDistanceTable();
    const Lattice& lat = arg_cluster->getLattice();
    // add copy of argument atom to calculate self contribution
    Atom_t arg_atom_uc = shiftToUnitCell(*arg_atom);
    cluster_atoms.back() = &arg_atom_uc;
    int lsqidx = 0;
    vector<double>::iterator fii = lsq_fi.begin();
    vector<double>::iterator dii = lsq_di.begin();
    vector<double>::iterator wtii = lsq_wt.begin();
    for (AtomSequenceIndex seq(cluster_atoms); !seq.finished(); seq.next())
    {
        R3::Vector rc_dd;
        double wt = sqrt(costToFitness( seq.ptr()->Badness() ));
        for (sph->rewind(); !sph->finished(); sph->next())
        {
            rc_dd = seq.ptr()->r + lat.cartesian(sph->mno) - arg_atom_uc.r;
            double d = R3::norm(rc_dd);
            if (d < rmin || d > rmax || d == 0.0)   continue;
            lsqidx++;
            assert(fii < lsq_fi.end());
            assert(dii < lsq_di.end());
            assert(wtii < lsq_wt.end());
            *(dii++) = d;
            *(wtii++) = wt;
            size_t nearidx = nearDistanceIndex(d);
            const double& dnear = dtgt[nearidx];
            *(fii++) = wt * (dnear - d);
            for (int j = 0; j < int(lsqParametersSize()); ++j)
            {
                _lsq_jacobian(lsqidx, j) = wt * (rc_dd[j]) / d;
            }
        }
    }
    cluster_atoms.back() = NULL;
    return lsq_fi;
}

double AtomCostCrystal::lsqJacobianGet(size_t m, size_t n) const
{
    // make sure all lsq_vectors are cached
    if (lsq_fi.empty())	    lsqComponents();
    return _lsq_jacobian(int(m), int(n));
}

// protected methods

void AtomCostCrystal::resizeArrays()
{
    partial_costs.assign(arg_cluster->NAtoms() + 1, 0.0);
}

Atom_t AtomCostCrystal::shiftToUnitCell(const Atom_t& a) const
{
    const Lattice& lat = arg_cluster->getLattice();
    // evaluate cartesian offset of argument atom
    R3::Vector fc_a = lat.fractional(a.r);
    // copy argument atom and shift it by floored fc_a
    Atom_t auc(a);
    auc.r = a.r - lat.cartesian(floor(fc_a));
    return auc;
}
