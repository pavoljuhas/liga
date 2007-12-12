/***********************************************************************
* Short Title: functions for obtaining subversion information
*
* Comments: subversion information definitions
*
* $Id: RegisterSVNId.hpp 1106 2007-07-04 07:12:12Z juhas $
* 
* <license text>
***********************************************************************/

#ifndef VERSION_HPP_INCLUDED
#define VERSION_HPP_INCLUDED

#include <string>

namespace Version {

using namespace std;

int getRevisionNumber();
const string& getRevision();
const string& getDate();
const string& getAuthor();
const string& getId();

}   // namespace Version

#endif	// VERSION_HPP_INCLUDED
