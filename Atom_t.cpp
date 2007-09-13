/***********************************************************************
* Short Title: definitions for Atom_t class
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include "Counter.hpp"
#include "Atom_t.hpp"

RegisterSVNId Atom_t_cpp_id("$Id$");

using namespace std;

////////////////////////////////////////////////////////////////////////
// class Atom_t
////////////////////////////////////////////////////////////////////////

// constructors

Atom_t::Atom_t(double rx0, double ry0, double rz0, double bad0) :
    fixed(false), ttp(LINEAR), badness(bad0)
{
    r[0] = rx0;
    r[1] = ry0;
    r[2] = rz0;
}

// public methods

double Atom_t::Badness() const
{
    return badness;
}

double Atom_t::FreeBadness() const
{
    return fixed ? 0.0 : badness;
}

double Atom_t::IncBadness(double db)
{
    badness += db;
    return badness;
}

double Atom_t::DecBadness(double db)
{
    badness -= db;
    if (badness < 0.0)	badness = 0.0;
    return badness;
}

double Atom_t::ResetBadness(double b)
{
    badness = b;
    return badness;
}

// non-member operators and functions

bool operator==(const Atom_t& a1, const Atom_t& a2)
{
    return equal(a1.r.data(), a1.r.data() + 3, a2.r.data());
}

double dist2(const Atom_t& a1, const Atom_t& a2)
{
    static Counter* distance_calls = Counter::getCounter("distance_calls");
    distance_calls->count();
    double dr[3] = { a1.r[0]-a2.r[0], a1.r[1]-a2.r[1], a1.r[2]-a2.r[2] };
    return dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
}

// End of file
