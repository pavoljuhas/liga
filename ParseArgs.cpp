/***********************************************************************
* Short Title: program arguments parser
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

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
}

void ParseArgs::List()
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
    while (1)
    {
	int this_option_optind = optind ? optind : 1;
	int option_index = 0;
	int c = getopt(argc, argv, optstring);
	if (c == -1)
	    break;
	else if (c == '?')
	    throw InvalidOption();
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
    while (1)
    {
	int this_option_optind = optind ? optind : 1;
	int option_index = 0;
	int c = getopt_long(argc, argv, optstring,
		longopts, &option_index);
	if (c == -1)
	    break;
	else if (c == '?')
	    throw InvalidOption();
	else
	{
	    const char *o = longopts[option_index].name;
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

