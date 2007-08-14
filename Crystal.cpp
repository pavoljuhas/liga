#include "/u24/juhas/programs/include/dbprint.h"
/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <sstream>
#include <boost/foreach.hpp>
#include "Crystal.hpp"
#include "Lattice.hpp"

using namespace std;

RegisterSVNId Crystal_cpp_id("$Id$");

////////////////////////////////////////////////////////////////////////
// class Crystal
////////////////////////////////////////////////////////////////////////

// constructors

Crystal::Crystal() : Molecule()
{ }

Crystal::Crystal(const DistanceTable& dtab) : Molecule(dtab)
{ }

Crystal::Crystal(const Crystal& crs) : Molecule(crs)
{ }

Crystal::~Crystal()
{ }

// public methods

Crystal& Crystal::operator=(const Crystal& crs)
{
    if (this == &crs)   return *this;
    // use base class copy
    Molecule& molthis = *this;
    molthis = crs;
    // copy Crystal specific data
    _lattice = crs._lattice;
    _rmin = crs._rmin;
    _rmax = crs._rmax;
    return *this;
}

void Crystal::setLattice(const Lattice& lat)
{
    _lattice.reset(new Lattice(lat));
}

const Lattice& Crystal::getLattice() const
{
    return *_lattice;
}

void Crystal::setRRange(double rmin, double rmax)
{
    _rmin = rmin;
    _rmax = rmax;
}

pair<double,double> Crystal::getRRange() const
{
    return make_pair(_rmin, _rmax);
}

// r-range extended by circum diameter of the primitive cell
pair<double,double> Crystal::getRExtent() const
{
    R3::Vector center(0.0, 0.0, 0.0);
    BOOST_FOREACH (Atom_t* ai, atoms)
    {
        center += ai->r;
    }
    center /= NAtoms();
    double max_offcenter = 0.0;
    BOOST_FOREACH (Atom_t* ai, atoms)
    {
        double ri = R3::norm(ai->r - center);
        if (ri > max_offcenter)     max_offcenter = ri;
    }
    // adjust for round-off errors, when not empty
    if (NAtoms() != 0)
    {
        const double epsd = sqrt(DOUBLE_EPS);
        max_offcenter = max_offcenter*(1.0 + epsd) + epsd;
    }
    double circum_diameter = 2*max_offcenter;
    return make_pair(_rmin - circum_diameter, _rmax + circum_diameter);
}

// private methods

void Crystal::init()
{
    DBPRINT("in Crystal::init()");
    _rmin = 0.0;
    _rmax = DOUBLE_MAX;
}

// End of file
