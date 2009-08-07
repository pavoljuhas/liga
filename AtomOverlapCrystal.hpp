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
* class AtomOverlap
*
* Comments: OverlapCost calculation for a Molecule
*
* $Id$
*
***********************************************************************/

#ifndef ATOMOVERLAPCRYSTAL_HPP_INCLUDED
#define ATOMOVERLAPCRYSTAL_HPP_INCLUDED

#include "AtomCostCrystal.hpp"

class AtomOverlapCrystal : public AtomCostCrystal
{
    public:

	// constructor
	AtomOverlapCrystal(const Crystal* crst);

	// public methods - overloaded
        virtual void resetFor(const Molecule* crst);
        virtual std::pair<double,int> pairCostCount(const R3::Vector& cv);

};  // class AtomOverlap

#endif  // ATOMOVERLAPCRYSTAL_HPP_INCLUDED
