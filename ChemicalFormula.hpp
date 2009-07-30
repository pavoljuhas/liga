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
* ChemicalFormula
*
* Comments: Storage and manipulation of ChemicalFormula.
*
* $Id$
*
***********************************************************************/

#ifndef CHEMICALFORMULA_HPP_INCLUDED
#define CHEMICALFORMULA_HPP_INCLUDED

#include <utility>
#include <list>
#include <string>

#include <boost/foreach.hpp>

typedef std::list< std::pair<std::string,int> > ChemicalFormula;

ChemicalFormula chemicalFormulaFromString(const std::string&);
std::string chemicalFormulaAsString(const ChemicalFormula&);

#endif  // CHEMICALFORMULA_HPP_INCLUDED
