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

	// class methods
	static std::string getFile(const std::string& id);
	static int getRevision(const std::string& id);
	static std::string getDate(const std::string& id);
	static std::string getAuthor(const std::string& id);
	static int lastRevision();
	static std::string lastDate();
	static std::string lastAuthor();
	static std::string lastId();

	// constructor
	RegisterSVNId(const char* id);

    private:

	// class data
	static int last_revision;
	static std::string last_date;
	static std::string last_author;

	// class methods
	static void processRecords();
	static std::map<std::string,std::string>& records(const char* id=NULL);

};  // class RegisterSVNId

#endif	// REGISTERSVNID_HPP_INCLUDED
