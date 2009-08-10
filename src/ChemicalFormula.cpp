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
#include <cassert>

#include "ChemicalFormula.hpp"

using namespace std;

// Constructors --------------------------------------------------------------

ChemicalFormula::ChemicalFormula()
{
}


ChemicalFormula::ChemicalFormula(const string& formula)
{
    this->fromString(formula);
}

// Public Methods ------------------------------------------------------------

int ChemicalFormula::countElements() const
{
    int cnt = 0;
    for (const_iterator sc = this->begin(); sc != this->end(); ++sc)
    {
        cnt += max(0, sc->second);
    }
    return cnt;
}


vector<string> ChemicalFormula::expand() const
{
    vector<string> rv;
    for (const_iterator sc = this->begin(); sc != this->end(); ++sc)
    {
        int c = max(0, sc->second);
        rv.insert(rv.end(), c, sc->first);
    }
    assert(int(rv.size()) == this->countElements());
    return rv;
}


void ChemicalFormula::fromString(const std::string& s)
{
    ChemicalFormula chfm;
    for (string::const_iterator c = s.begin(); c != s.end(); ++c)
    {
        if (isblank(*c))  continue;
        if (chfm.empty() || isupper(*c))
        {
            chfm.push_back(ChemicalFormula::value_type("", 0));
        }
        string& smblcnt = chfm.back().first;
        smblcnt += *c;
    }
    for (iterator sc = chfm.begin(); sc != chfm.end(); ++sc)
    {
        string& smblcnt = sc->first;
        string::size_type pcnt;
        pcnt = smblcnt.find_last_not_of("0123456789") + 1;
        string scnt = smblcnt.substr(pcnt);
        smblcnt.erase(pcnt);
        sc->second = scnt.empty() ? 1 : atoi(scnt.c_str());
    }
    // assign when everything went 
    this->swap(chfm);
}


string ChemicalFormula::toString(string separator) const
{
    ostringstream rv;
    for (const_iterator sc = this->begin(); sc != this->end(); ++sc)
    {
        if (sc != this->begin())  rv << separator;
        if (sc->second > 0)  rv << sc->first;
        if (sc->second > 1)  rv << sc->second;
    }
    return rv.str();
}

// End of file
