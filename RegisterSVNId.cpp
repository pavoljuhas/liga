/***********************************************************************
* Short Title: class for managing SVN Id strings from each object
*
* Comments: Useful for obtaining executable version.  Cannot register
*	    itself, but this is a small deal.
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <sstream>
#include "RegisterSVNId.hpp"

using namespace std;


////////////////////////////////////////////////////////////////////////
// class RegisterSVNId
////////////////////////////////////////////////////////////////////////

// class data

map<string,string> RegisterSVNId::svnids;
int RegisterSVNId::last_revision = 0;
string RegisterSVNId::last_date;
string RegisterSVNId::last_author;
string RegisterSVNId::last_id;

// public class methods

string RegisterSVNId::getFile(const string& id)
{
    istringstream idline(id);
    string ignore, wfile;
    idline >> ignore >> wfile;
    return wfile;
}

string RegisterSVNId::getDate(const string& id)
{
    istringstream idline(id);
    string ignore, wdate, wtime;
    idline >> ignore >> ignore >> ignore >> wdate >> wtime;
    return wdate + " " + wtime;
}

int RegisterSVNId::getRevision(const string& id)
{
    istringstream idline(id);
    string ignore;
    int rev = 0;
    idline >> ignore >> ignore >> rev;
    return rev;
}

string RegisterSVNId::getAuthor(const string& id)
{
    istringstream idline(id);
    string ignore, wauthor;
    idline >> ignore >> ignore >> ignore >> ignore >> ignore >> wauthor;
    return wauthor;
}

// constructor

RegisterSVNId::RegisterSVNId(const char* id)
{
    /*
    string myid = id;
    svnids[getFile(myid)] = myid;
    if (getRevision(myid) > last_revision)
    {
	last_revision = getRevision(myid);
	last_date = getDate(myid);
	last_author = getAuthor(myid);
	ostringstream idstrm;
	idstrm << last_revision << ' ' << last_date << ' ' << last_author;
	last_id = idstrm.str();
    }
    */
}

// End of RegisterSVNId.cpp
