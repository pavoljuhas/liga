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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <getopt.h>
#include "ioerror.hpp"

using namespace std;

struct ParseArgsError : public runtime_error
{
    ParseArgsError (string what_arg = "") :
	runtime_error(what_arg) { }
};

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
    void Dump();
    void ReadPars(const char *file);
    istream& ReadPars(istream& fid = cin);
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
	ostringstream oss;
	oss << "parameter '" << par << "' is not defined";
	throw ParseArgsError(oss.str());
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
	ostringstream oss;
	oss << "parameter '" << par << "' is not defined";
	throw ParseArgsError(oss.str());
    }
    vector<T> v;
    istringstream iss(pars[par]);
    T val;
    while (iss >> val)
	v.push_back(val);
    return v;
}

#endif
