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
    template <class T> T GetPar(string par);
    template <class T> T GetPar(string par, T defval);
    template <typename T> vector<T> GetParVec(string par);
private:
    void init();
    void do_getopt();
    void do_getopt_long();
    void arg_or_par(const char *s);
};

template<class T> T ParseArgs::GetPar(string par)
{
    if (!pars.count(par))
    {
	cerr << "parameter '" << par << "' is not defined" << endl;
	throw InvalidParameter();
    }
    T val;
    istringstream iss(pars[par]);
    iss >> val;
    return val;
}

template<class T> T ParseArgs::GetPar(string par, T defval)
{
    return pars.count(par) ? GetPar<T>(par) : defval;
}

template <typename T> vector<T> GetParVec(string par)
{
    if (!pars.count(par))
    {
	cerr << "parameter '" << par << "' is not defined" << endl;
	throw InvalidParameter();
    }
    vector<T> v;
    istringstream iss(pars[par]);
    T val;
    while (iss >> val)
	v.push_back(val);
    return v;
}

#endif
