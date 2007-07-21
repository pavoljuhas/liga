/***********************************************************************
* Short Title: derived class from Molecule for periodic crystal structure
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef CRYSTAL_T_HPP_INCLUDED
#define CRYSTAL_T_HPP_INCLUDED

#include "BGAlib.hpp"
#include "RegisterSVNId.hpp"

class Lattice;

namespace {
RegisterSVNId Crystal_hpp_id("$Id$");
}

class Crystal : public Molecule
{
    public:

	// data

	// constructors
	Crystal();
	Crystal(const DistanceTable&);
	Crystal(const Crystal&);
	Crystal& operator=(const Crystal&);
	~Crystal();

        // methods

    private:

	// crystal specific data
	Lattice* lattice;
	double rmin;
	double rmax;

        // methods
        void init();



};


#endif	// CRYSTAL_T_HPP_INCLUDED
