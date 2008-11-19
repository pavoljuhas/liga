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

#include <string>
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
        friend class Crystal;

        // constructors
	Atom_t(const std::string& elsmbl,
                double rx0, double ry0, double rz0, double bad0=0.0);
        template <class V>
            Atom_t(const std::string& elsmbl, const V& r0, double bad0=0.0);

        // data
        const std::string element;
        R3::Vector r;
	bool fixed;
	triangulation_type ttp;

        // methods
	const double& Badness() const;
	double FreeBadness() const;
	void IncBadness(const double& db);
	void DecBadness(const double& db);
	void ResetBadness(double b=0.0);

    private:

        // data
	double _badness;
	mutable int pmxidx;     // pair matrix index
};

// declarations of non-member operators and functions

bool operator==(const Atom_t& a1, const Atom_t& a2);

////////////////////////////////////////////////////////////////////////
// template and inline functions
////////////////////////////////////////////////////////////////////////

// template constructor

template <class V>
Atom_t::Atom_t(const std::string& elsmbl, const V& r0, double bad0) :
    element(elsmbl), fixed(false), ttp(LINEAR), _badness(bad0)
{
    r = r0[0], r0[1], r0[2];
}


#endif	// ATOM_T_HPP_INCLUDED
