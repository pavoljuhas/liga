/***********************************************************************
* Short Title: class for registering SVN Id strings from each object
*
* Comments: Useful for obtaining executable version.  Cannot register
*	    itself, but this is a small deal.
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef REGISTERSVNID_HPP_INCLUDED
#define REGISTERSVNID_HPP_INCLUDED

#include <string>
#include <map>

class RegisterSVNId
{
    public:

	// class data
	static std::map<std::string,std::string> svnids;
	static int last_revision;
	static std::string last_date;
	static std::string last_author;
	static std::string last_id;

	// class methods
	static std::string getFile(const std::string& id);
	static int getRevision(const std::string& id);
	static std::string getDate(const std::string& id);
	static std::string getAuthor(const std::string& id);

	// constructor
	RegisterSVNId(const char* idstring);

    private:

	// data
	const std::string myid;

};  // class RegisterSVNId

#endif	// REGISTERSVNID_HPP_INCLUDED
