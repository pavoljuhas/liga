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
    cmd_h = argv[0];
    string::size_type pslash = cmd_h.find_last_of('/');
    cmd_t = cmd_h.substr(pslash+1);
    cmd_h.erase(pslash);
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

void ParseArgs::do_getopt()
{
    while (1)
    {
	int this_option_optind = optind ? optind : 1;
	int option_index = 0;
	int c = getopt(argc, argv, optstring);
	if (c == -1)
	    break;
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
	const char *o = longopts[option_index].name;
	opts[o] = (optarg) ? optarg : "";
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
	pars[string(s, peq-1-s)] = string(peq+1);
    else
	args.push_back(s);
}

void ParseArgs::ReadPars(const char *file)
{}
