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
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <getopt.h>
#include "BGAutils.hpp"

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
    template<typename T> T GetPar(string par);
    std::vector<int> ExpandRangePar(string par);
    template<typename T> T GetPar(string par, T defval);
    template<typename T> vector<T> GetParVec(string par);
    void ValidatePars(list<string>& validpars);
    inline bool ispar(const string& p) { return pars.count(p); }
    inline bool isopt(const string& o) { return opts.count(o); }
private:
    void init();
    void do_getopt();
    void do_getopt_long();
    void arg_or_par(const char *s);
    map<string,bool> cmdl_par;
};

template<typename T> T ParseArgs::GetPar(string par)
{
    if (!pars.count(par))
    {
	ostringstream oss;
	oss << "parameter '" << par << "' is not defined";
	throw ParseArgsError(oss.str());
    }
    T val;
    istringstream iss(pars[par]);
    if ( !(iss >> val) )
    {
	ostringstream oss;
	oss << "invalid value for parameter '" << par << "'";
	throw ParseArgsError(oss.str());
    }
    return val;
}

template<typename T> T ParseArgs::GetPar(string par, T defval)
{
    return pars.count(par) ? GetPar<T>(par) : defval;
}

template<typename T> vector<T> ParseArgs::GetParVec(string par)
{
    if (!pars.count(par))
    {
	ostringstream oss;
	oss << "parameter '" << par << "' is not defined";
	throw ParseArgsError(oss.str());
    }
    // replace all commas in par with <space>
    string values(pars[par]);
    for (   string::size_type pcomma = values.find(',');
	    pcomma != string::npos; pcomma = values.find(',', pcomma) )
    {
	values[pcomma] = ' ';
    }
    vector<T> v;
    istringstream iss(values);
    T val;
    while (iss >> val)
	v.push_back(val);
    return v;
}

#endif
