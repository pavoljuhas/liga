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

using namespace std;

////////////////////////////////////////////////////////////////////////
// class Atom_t
////////////////////////////////////////////////////////////////////////

// constructors

Atom_t::Atom_t(double rx0, double ry0, double rz0, double bad0) :
    fixed(false), ttp(LINEAR), _badness(bad0)
{
    r[0] = rx0;
    r[1] = ry0;
    r[2] = rz0;
}

// public methods

const double& Atom_t::Badness() const
{
    return this->_badness;
}

double Atom_t::FreeBadness() const
{
    double rv = fixed ? 0.0 : this->Badness();
    return rv;
}

void Atom_t::IncBadness(const double& db)
{
    this->_badness += db;
}

void Atom_t::DecBadness(const double& db)
{
    this->_badness -= db;
    if (this->_badness < 0.0)	this->ResetBadness();
}

void Atom_t::ResetBadness(double b)
{
    this->_badness = b;
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
