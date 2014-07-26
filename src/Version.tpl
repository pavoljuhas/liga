/***********************************************************************
* Short Title: functions for obtaining subversion information
*
* Comments: subversion information definitions.
*     DO NOT edit Version.cpp, change Version.tpl instead.
*
* <license text>
***********************************************************************/

#include <sstream>
#include "Version.hpp"

namespace NS_VERSION {

int getRevisionNumber()
{
    static int rev = 0;
    if (!rev)   istringstream (getRevision()) >> rev;
    return rev;
}

const string& getRevision()
{
    static const string srev = "%(LastChangedRev)s";
    return srev;
}

const string& getDate()
{
    static const string sdate = "%(LastChangedDate)s";
    return sdate;
}

const string& getAuthor()
{
    static const string sauthor = "%(LastChangedAuthor)s";
    return sauthor;
}

const string& getId()
{
    static string sid;
    if (sid.empty())
    {
        ostringstream idstrm;
        idstrm << getRevision() << ' ' << getDate() << ' ' << getAuthor();
        sid = idstrm.str();
    }
    return sid;
}

}   // namespace NS_VERSION

// vim:ft=cpp:
// End of Version.cpp
