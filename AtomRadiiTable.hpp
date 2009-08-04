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
* AtomRadiiTable
*
* Comments: Storage of empirical atomic radii
*
* $Id$
*
***********************************************************************/

#ifndef ATOMRADIITABLE_HPP_INCLUDED
#define ATOMRADIITABLE_HPP_INCLUDED

#include <map>
#include <string>
#include <sstream>
#include <stdexcept>


class AtomRadiiTable : public std::map<std::string,double>
{
    public:

        /// fast value lookup, which does not change the table.
        double lookup(const std::string& smbl) const
        {
            using namespace std;
            const_iterator ai = this->find(smbl);
            if (ai == this->end())
            {
                ostringstream emsg;
                emsg << "Undefined radius for element '" << smbl << "'.";
                throw invalid_argument(emsg.str());
            }
            return ai->second;
        }
};


#endif  // ATOMRADIITABLE_HPP_INCLUDED
