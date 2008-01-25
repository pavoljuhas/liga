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
#include "Atom_t.hpp"
#include "AtomCostCrystal.hpp"
#include "LigaUtils.hpp"
#include "PointsInSphere.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// class Crystal
////////////////////////////////////////////////////////////////////////


// constructors

Crystal::Crystal() : Molecule()
{
    init();
}

Crystal::Crystal(const DistanceTable& dtab) : Molecule(dtab)
{
    init();
}

Crystal::Crystal(const Crystal& crs) : Molecule(crs)
{
    init();
}

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
    this->_lattice = crs._lattice;
    this->_rmin = crs._rmin;
    this->_rmax = crs._rmax;
    this->_cost = crs._cost;
    this->_cost_cached = crs._cost_cached;
    this->_cost_of_lattice = crs._cost_of_lattice;
    this->_cost_of_lattice_cached = crs._cost_of_lattice_cached;
    return *this;
}

void Crystal::setDistanceTable(const DistanceTable& dtbl)
{
    vector<double> dtblunique = dtbl.unique();
    this->_distance_table.reset(new DistanceTable(dtblunique));
    this->uncacheCost();
}

const DistanceTable& Crystal::getDistanceTable() const
{
    return *_distance_table;
}

void Crystal::setLattice(const Lattice& lat)
{
    _lattice.reset(new Lattice(lat));
    this->uncacheCost();
}

const Lattice& Crystal::getLattice() const
{
    return *_lattice;
}

void Crystal::setRRange(double rmin, double rmax)
{
    _rmin = rmin;
    _rmax = rmax;
    this->uncacheCost();
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

double Crystal::cost() const
{
    if (!_cost_cached)
    {
        Recalculate();
        _cost = NormBadness() + costOfLattice();
        _cost_cached = true;
    }
    return _cost;
}

double Crystal::costOfLattice() const
{
    if (!_cost_of_lattice_cached)
    {
        _cost_of_lattice = 0.0;
        if (_lattice.get())
        {
            AtomCostCrystal* atomcost = static_cast<AtomCostCrystal*>(
                    this->getAtomCostCalculator());
            _cost_of_lattice = atomcost->costOfLattice();
        }
        _cost_of_lattice_cached = true;
    }
    return _cost_of_lattice;
}

// protected methods

AtomCost* Crystal::getAtomCostCalculator() const
{
    static AtomCostCrystal the_acc(this);
    the_acc.resetFor(this);
    return &the_acc;
}

// private methods

void Crystal::init()
{
    this->_rmin = 0.0;
    this->_rmax = 0.0;
    this->uncacheCost();
}


void Crystal::uncacheCost()
{
    this->_cost_cached = false;
    this->_cost_of_lattice_cached = false;
}


// End of file
