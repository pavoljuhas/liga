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

#include <cstdlib>
#include <cctype>
#include <sstream>
#include <boost/foreach.hpp>
#include "ChemicalFormula.hpp"

using namespace std;


ChemicalFormula chemicalFormulaFromString(const std::string& sfml)
{
    ChemicalFormula rv;
    BOOST_FOREACH (string::value_type c, sfml)
    {
        if (isblank(c))  continue;
        if (rv.empty() || isupper(c))
        {
            rv.push_back(ChemicalFormula::value_type("", 0));
        }
        string& smblcnt = rv.back().first;
        smblcnt += c;
    }
    BOOST_FOREACH (ChemicalFormula::value_type& fmitem, rv)
    {
        string& smblcnt = fmitem.first;
        string::size_type pcnt;
        pcnt = smblcnt.find_last_not_of("0123456789") + 1;
        string scnt = smblcnt.substr(pcnt);
        smblcnt.erase(pcnt);
        fmitem.second = scnt.empty() ? 1 : atoi(scnt.c_str());
    }
    return rv;
}


string chemicalFormulaAsString(const ChemicalFormula& formula)
{
    ostringstream rv;
    BOOST_FOREACH (const ChemicalFormula::value_type& fmitem, formula)
    {
        rv << fmitem.first;
        if (fmitem.second != 1)   rv << fmitem.second;
    }
    return rv.str();
}

// End of file
