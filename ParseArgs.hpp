/***********************************************************************
* Short Title: program arguments parser
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef PARSEARGS_HPP_INCLUDED
#define PARSEARGS_HPP_INCLUDED
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <getopt.h>

// exceptions
struct IOError { };
struct InvalidOption { };
struct InvalidParameter { };

using namespace std;

class ParseArgs
{
public:
    // constructor
    ParseArgs(int nargc, char * const nargv[]);
    ParseArgs(int nargc, char * const nargv[], char *optstring);
    ParseArgs(int nargc, char * const nargv[], char *optstring,
	    const struct option *longopts);
    int argc;
    char * const * argv;
    char *optstring;
    const struct option *longopts;
    map<string,string> opts;
    map<string,string> pars;
    vector<string> args;
    string cmd_h, cmd_t;
    void Parse();
    void List();
    void ReadPars(const char *file);
    template <class T> T getpar(string p);
    template <class T> T getpar(string p, T v);
private:
    void init();
    void do_getopt();
    void do_getopt_long();
    void arg_or_par(const char *s);
};

template<class T> T ParseArgs::getpar(string p)
{
    if (!pars.count(p))
    {
	cerr << "parameter '" << p << "' is not defined" << endl;
	throw InvalidParameter();
    }
    T v;
    istringstream iss(pars[p]);
    iss >> v;
    return v;
}

template<class T> T ParseArgs::getpar(string p, T v)
{
    return pars.count(p) ? getpar<T>(p) : v;
}

#endif
