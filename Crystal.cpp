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

Crystal& Crystal::operator=(const Crystal& crs)
{
    if (this == &crs)   return *this;
    // use base class copy
    Molecule& molthis = *this;
    molthis = crs;
    // copy Crystal specific data
    lattice = crs.lattice;
    rmin = crs.rmin;
    rmax = crs.rmax;
    return *this;
}

Crystal::~Crystal()
{ }


// private methods

void Crystal::init()
{
    DBPRINT("in Crystal::init()");
    rmin = 0.0;
    rmax = DOUBLE_MAX;
    lattice = NULL;
}

// End of file
