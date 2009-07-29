/***********************************************************************
* Short Title: definitions for Atom_t class
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include "LigaUtils.hpp"
#include "Atom_t.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// class Atom_t
////////////////////////////////////////////////////////////////////////

// constructors

Atom_t::Atom_t(const string& elsmbl,
        double rx, double ry, double rz) : element(elsmbl)
{
    this->init(rx, ry, rz);
}

// public methods

const double& Atom_t::Badness() const
{
    return this->_badness;
}

double Atom_t::FreeBadness() const
{
    double rv = this->fixed ? 0.0 : this->Badness();
    return rv;
}

void Atom_t::IncBadness(const double& db)
{
    this->_badness += db;
}

void Atom_t::DecBadness(const double& db)
{
    this->_badness -= db;
    // take care of round-off errors
    if (this->_badness < NS_LIGA::eps_cost)     this->ResetBadness(0.0);
}

void Atom_t::ResetBadness(double b)
{
    this->_badness = b;
}

// private methods

void Atom_t::init(double rx, double ry, double rz)
{
    this->r = rx, ry, rz;
    this->fixed = false;
    this->ttp = LINEAR;
    this->ResetBadness(0.0);
}

// non-member operators and functions

bool operator==(const Atom_t& a1, const Atom_t& a2)
{
    bool rv = (a1.element == a2.element &&
            equal(a1.r.data(), a1.r.data() + 3, a2.r.data()));
    return rv;
}

// End of file
