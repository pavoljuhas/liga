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
#include "RegisterSVNId.hpp"

namespace {
RegisterSVNId ParseArgs_hpp_id("$Id$");
}

class ParseArgsError : public std::runtime_error
{
    public:
	ParseArgsError (std::string what_arg="") :
	    std::runtime_error(what_arg)
	{ }
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
    std::map<std::string,std::string> opts;
    std::map<std::string,std::string> pars;
    std::vector<std::string> args;
    std::string cmd_h, cmd_t;
    void Parse();
    void Dump();
    void ReadPars(const char *file);
    std::istream& ReadPars(std::istream& fid=std::cin);
    template<typename T> T GetPar(std::string par);
    std::vector<int> ExpandRangePar(std::string par);
    template<typename T> T GetPar(std::string par, T defval);
    template<typename T> std::vector<T> GetParVec(std::string par);
    void ValidatePars(const std::list<std::string>& validpars);
    inline bool ispar(const std::string& p) { return pars.count(p); }
    inline bool isopt(const std::string& o) { return opts.count(o); }
private:
    void init();
    void do_getopt();
    void do_getopt_long();
    void arg_or_par(const char *s);
    std::map<std::string,bool> cmdl_par;
};

template<typename T> T ParseArgs::GetPar(std::string par)
{
    if (!pars.count(par))
    {
	std::ostringstream oss;
	oss << "parameter '" << par << "' is not defined";
	throw ParseArgsError(oss.str());
    }
    T val;
    std::istringstream iss(pars[par]);
    if ( !(iss >> val) )
    {
	std::ostringstream oss;
	oss << "invalid value for parameter '" << par << "'";
	throw ParseArgsError(oss.str());
    }
    return val;
}

template<typename T> T ParseArgs::GetPar(std::string par, T defval)
{
    return pars.count(par) ? GetPar<T>(par) : defval;
}

template<typename T> std::vector<T> ParseArgs::GetParVec(std::string par)
{
    if (!pars.count(par))
    {
	std::ostringstream oss;
	oss << "parameter '" << par << "' is not defined";
	throw ParseArgsError(oss.str());
    }
    // replace all commas in par with <space>
    std::string values(pars[par]);
    for (   std::string::size_type pcomma = values.find(',');
	    pcomma != std::string::npos; pcomma = values.find(',', pcomma) )
    {
	values[pcomma] = ' ';
    }
    std::vector<T> v;
    std::istringstream iss(values);
    T val;
    while (iss >> val)
	v.push_back(val);
    return v;
}

#endif
