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
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <getopt.h>

// exceptions
struct IOError { };
struct InvalidOption { };

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
private:
    void init();
    void do_getopt();
    void do_getopt_long();
    void arg_or_par(const char *s);
};

#endif
