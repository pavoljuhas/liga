/***********************************************************************
* Short Title: program arguments parser
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <iostream>
#include "ParseArgs.hpp"

ParseArgs::ParseArgs(int nargc, char * const nargv[]) :
    argc(nargc), argv(nargv)
{
    optstring = NULL;
    longopts = NULL;
    init();
}

ParseArgs::ParseArgs(int nargc, char * const nargv[], char *noptstring) :
    argc(nargc), argv(nargv), optstring(noptstring)
{
    longopts = NULL;
    init();
}

ParseArgs::ParseArgs(int nargc, char * const nargv[], char *noptstring,
	    const struct option *nlongopts) :
    argc(nargc), argv(nargv), optstring(noptstring), longopts(nlongopts)
{
    init();
}

void ParseArgs::Parse()
{
    if (longopts && optstring)
	do_getopt_long();
    else if (optstring)
	do_getopt();
    else 
    {
	for (int i = 1; i < argc; ++i)
	    arg_or_par(argv[i]);
    }
}

void ParseArgs::ReadPars(const char *file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	ostringstream oss;
	oss << "Unable to read '" << file << "'";
	throw IOError(oss.str());
    }
    ReadPars(fid);
    fid.close();
}

istream& ParseArgs::ReadPars(istream& fid)
{
    const char *blank = " \t\n";
    string line;
    for (int nr = 1; !fid.eof() && getline(fid, line); ++nr)
    {
	string::size_type lb, le;
	lb = line.find_first_not_of(blank);
	le = line.find_last_not_of(blank)+1;
	line = (lb == string::npos) ? "" : line.substr(lb, le-lb);
	if (line.length() == 0 || line[0] == '#')
	    continue;
	string::size_type eq, pe, vb;
	eq = line.find('=');
	if (eq == string::npos)
	{
	    ostringstream oss;
	    oss << nr << ": missing equal symbol";
	    throw ParseArgsError(oss.str());
	}
	pe = line.find_last_not_of(string(blank) + "=", eq);
	pe = (pe == string::npos) ? eq : pe+1;
	string par = line.substr(0, pe);
	bool ispar = isalpha(par[0]);
	for (   string::iterator ii = par.begin();
		ii != par.begin() && ispar; ++ii)
	{
	    ispar = isalnum(*ii);
	}
	if (!ispar)
	{
	    ostringstream oss;
	    oss << nr << ": invalid parameter name '" << par << "'";
	    throw ParseArgsError(oss.str());
	}
	vb = line.find_first_not_of(blank, eq+1);
	pars[par] = (vb == string::npos) ? "" : line.substr(vb);
    }
    return fid;
}

void ParseArgs::Dump()
{
    cout << "cmd_h = '" << cmd_h << "'" << endl;
    cout << "cmd_t = '" << cmd_t << "'" << endl;
    typedef map<string,string>::iterator MSSit;
    for (MSSit ii = opts.begin(); ii != opts.end(); ++ii)
	cout << "opts[" << ii->first << "] = '" << ii->second << "'" << endl;
    for (MSSit ii = pars.begin(); ii != pars.end(); ++ii)
	cout << "pars[" << ii->first << "] = '" << ii->second << "'" << endl;
    for (int i = 0; i < args.size(); ++i)
	cout << "args[" << i << "] = '" << args[i] << "'" << endl;
}

void ParseArgs::init()
{
    cmd_h = argv[0];
    string::size_type pslash = cmd_h.find_last_of('/');
    if (pslash != string::npos)
    {
	cmd_t = cmd_h.substr(pslash+1);
	cmd_h.erase(pslash);
    }
    else
	cmd_h.swap(cmd_t);
}

void ParseArgs::do_getopt()
{
    while (true)
    {
	int this_option_optind = optind ? optind : 1;
	int option_index = 0;
	int c = getopt(argc, argv, optstring);
	if (c == -1)
	    break;
	else if (c == '?')
	    throw ParseArgsError();
	else
	    opts[string(1,c)] = (optarg) ? optarg : "";
    }
    for (; optind < argc; ++optind)
    {
	arg_or_par(argv[optind]);
    }
}

void ParseArgs::do_getopt_long()
{
    while (true)
    {
	int this_option_optind = optind ? optind : 1;
	int option_index = 0;
	int c = getopt_long(argc, argv, optstring,
		longopts, &option_index);
	if (c == -1)
	    break;
	string o;
	switch (c)
	{
	    case '?':
		throw ParseArgsError();
		break;
	    case 0:
		o = longopts[option_index].name;
		opts[o] = (optarg) ? optarg : "";
		break;
	    default:
		o = c;
		opts[o] = (optarg) ? optarg : "";
	}
    }
    for (; optind < argc; ++optind)
    {
	arg_or_par(argv[optind]);
    }
}

void ParseArgs::arg_or_par(const char *s)
{
    char *peq = strchr(s, '=');
    if (peq == NULL)
    {
	args.push_back(s);
	return;
    }
    // here it looks like a parameter
    bool ispar = isalpha(s[0]);
    for (const char *p = s+1; ispar && p < peq; ++p)
	ispar = isalnum(*p);
    if (ispar)
	pars[string(s, peq-s)] = string(peq+1);
    else
	args.push_back(s);
}
