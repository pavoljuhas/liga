#include "/u24/juhas/programs/include/dbprint.h"
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

int RegisterSVNId::lastRevision()
{
    processRecords();
    return last_revision;
}

string RegisterSVNId::lastDate()
{
    processRecords();
    return last_date;
}

string RegisterSVNId::lastAuthor()
{
    processRecords();
    return last_author;
}

string RegisterSVNId::lastId()
{
    processRecords();
    ostringstream idstrm;
    idstrm << last_revision << ' ' << last_date << ' ' << last_author;
    return idstrm.str();
}

// constructor

RegisterSVNId::RegisterSVNId(const char* id)
{
    records(id);
}

// private class data

int RegisterSVNId::last_revision = 0;
string RegisterSVNId::last_date;
string RegisterSVNId::last_author;

// private class methods

void RegisterSVNId::processRecords()
{
    static size_t last_records_size = 0;
    if (records().size() == last_records_size)	return;
    map<string,string>& svnids = records();
    last_records_size = svnids.size();
    map<string,string>::iterator ii;
    for (ii = svnids.begin(); ii != svnids.end(); ++ii)
    {
	if (getRevision(ii->second) > last_revision)
	{
	    last_revision = getRevision(ii->second);
	    last_date = getDate(ii->second);
	    last_author = getAuthor(ii->second);
	}
    }
}

// helper function for holding data during initialization
map<string,string>& RegisterSVNId::records(const char* id)
{
    static map<string,string> svnids;
    if (id)
    {
	string idstr = id;
	svnids[getFile(idstr)] = idstr;
    }
    return svnids;
}

// End of RegisterSVNId.cpp
