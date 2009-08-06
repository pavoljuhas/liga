/***********************************************************************
*
* Liga Algorithm    for structure determination from pair distances
*                   Pavol Juhas
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* class AtomOverlapCrystal
*
* Comments: OverlapCost calculation for a Molecule
*
* $Id$
*
***********************************************************************/

#include <vector>

#include "AtomOverlapCrystal.hpp"
#include "AtomSequence.hpp"
#include "Crystal.hpp"
#include "Lattice.hpp"
#include "LigaUtils.hpp"
#include "Counter.hpp"

using namespace std;

// Constructor ---------------------------------------------------------------

AtomOverlapCrystal::AtomOverlapCrystal(const Crystal* crst) :
    AtomCostCrystal(crst)
{ }

// Public Methods ------------------------------------------------------------

void AtomOverlapCrystal::resetFor(const Molecule* clust)
{
    this->AtomCost::resetFor(clust);
    this->arg_cluster = static_cast<const Crystal*>(clust);
    assert(this->arg_cluster == this->AtomCost::arg_cluster);
    assert(!this->use_distances);
    this->_rmax = 2 * arg_cluster->getMaxAtomRadius();
    pair<double,double> rext = arg_cluster->getRExtent();
    _sph.reset(new PointsInSphere(rext.first, rext.second,
                arg_cluster->getLattice()) );
}


pair<double,int>
AtomOverlapCrystal::pairCostCount(const R3::Vector& cv, bool skipzero)
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
        double r0r1 = arg_atom->radius + crst_atom->radius;
	double dd = (d < r0r1) ? (r0r1 - d) : 0.0;
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
    pair<double,int> rv(paircost, paircount);
    return rv;
}


// End of file
