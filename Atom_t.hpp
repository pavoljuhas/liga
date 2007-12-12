/***********************************************************************
* Short Title: object declarations for Atom_t class
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef ATOM_T_HPP_INCLUDED
#define ATOM_T_HPP_INCLUDED

#include "R3linalg.hpp"

// may be later moved to dedicated triangulation class
enum triangulation_type { LINEAR, PLANAR, SPATIAL, NTGTYPES };

////////////////////////////////////////////////////////////////////////
// class Atom_t
////////////////////////////////////////////////////////////////////////

class Atom_t
{
    public:

        // friends
        friend class Molecule;

        // constructors
	Atom_t(double rx0, double ry0, double rz0, double bad0=0.0);
        template <class V> Atom_t(const V& r0, double bad0=0.0);

        // data
        R3::Vector r;
	bool fixed;
	triangulation_type ttp;

        // methods
	double Badness() const;
	double FreeBadness() const;
	double IncBadness(double db);
	double DecBadness(double db);
	double ResetBadness(double b=0.0);

    private:

        // data
	double badness;
	int pmxidx;	    // pair matrix index
};

// declarations of non-member operators and functions

double dist(const Atom_t& a1, const Atom_t& a2);
bool operator==(const Atom_t& a1, const Atom_t& a2);
double dist2(const Atom_t& a1, const Atom_t& a2);

////////////////////////////////////////////////////////////////////////
// template and inline functions
////////////////////////////////////////////////////////////////////////

// template constructor

template <class V>
Atom_t::Atom_t(const V& r0, double bad0) :
    fixed(false), ttp(LINEAR), badness(bad0)
{
    r = r0[0], r0[1], r0[2];
}

// inline non-member functions

inline double dist(const Atom_t& a1, const Atom_t& a2)
{
    return sqrt(1.0*dist2(a1, a2));
}

#endif	// ATOM_T_HPP_INCLUDED
