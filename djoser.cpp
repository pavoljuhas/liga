/***********************************************************************
* Short Title: molecule reconstruction from distance table
*
* Comments: not particularly inteligent search for molecular configuration
*
* $Id$
***********************************************************************/

#include <limits>
#include <unistd.h>
#include "ioerror.hpp"
#include "ParseArgs.hpp"
#include "BGAlib.hpp"

void print_version()
{
    using namespace std;
    const char *ver =
	"$Id$"
#   if defined(__DATE__) && defined(__TIME__)
	"\ncompiled " __DATE__ " " __TIME__
#   endif
	"\n";
    cout << ver;
}

void print_help(ParseArgs& a)
{
	/*

usage: << a.cmd_t << [-p PAR_FILE] [DISTHIST] [par1=val1 par2=val2...]
run drunkwalk simulation using DISTHIST data.  Parameters can be set
in PAR_FILE or on the command line, which overrides PAR_FILE.
Options:
  -p, --parfile=FILE    read parameters from FILE
  -h, --help            display this message
  -v, --version         show program version
IO parameters:
  disthist=FILE         target distance table
  outstru=FILE          where to save the best full molecule
  inistru=FILE          initial structure [empty box]
  snapshot=FILE         live molecule structure
  snaprate=int          number of iterations between snapshot updates
  frames=FILE           save intermediate structures to FILE.iteration
  framesrate=int        number of iterations between frame saves
Walk parameters
  tol_dd=double         [inf] distance is not used when dd=|d-d0|>=tol_dd
  tol_bad=double        target value of full molecule badness
  logsize=int           [10] last steps used for success rate evaluation
  eprob_max=double      high limit of evolve probability
  eprob_min=double      low limit of evolve probability
  bustprob=double       probability of forcing the full structure built
  evolve_jump=bool      [true] allow additions of several atoms
  evolve_frac=double    selection badness threshold of tested atoms
  penalty=string        dd penalty function [pow2], fabs, well
  dist_trials=int       [10] good distance atoms to try
  tri_trials=int        [20] godd triangle atoms to try
  pyr_trials=int        [1000] good pyramid atoms to try


*/

}

int main(int argc, char *argv[])
{
    using namespace std;

    char *short_options =
	"p:hv";
    // parameters and options
    option long_options[] = {
	{"parfile", 1, 0, 'p'},
	{"help", 0, 0, 'h'},
	{"version", 0, 0, 'v'},
	{0, 0, 0, 0}
    };
    ParseArgs a(argc, argv, short_options, long_options);
    try {
	a.Parse();
    }
    catch (ParseArgsError) {
	return EXIT_FAILURE;
    }
    a.Dump();
    if (a.opts.count("h") || argc == 1)
    {
	print_help(a);
	return EXIT_SUCCESS;
    }
    else if (a.opts.count("v"))
    {
	print_version();
	return EXIT_SUCCESS;
    }





    const int logsize = 10;
    const double pemin = 0.25;
    const double pemax = 0.75;
    const double pallway = 0.01;
    const double avgmb = 0.01*0.01;
    const double tol_dd = 0.1;
    ////////////////////////////////////////////////////////////////////////
    if (argc == 1)
    {
	cerr << "usage: " <<
	    "molTest08 distance_file.dss [seed] [snapshot_file.xyz]" << endl;
	return EXIT_SUCCESS;
    }
    // here argc > 1
    char *distance_file = argv[1];
    DistanceTable *dtab;
    try
    {
	dtab = new DistanceTable(distance_file);
    }
    catch (IOError)
    {
	return EXIT_FAILURE;
    }
    char *snapshot_file = NULL;
    if (argc > 2 && strlen(argv[2]) > 0)
    {
	unsigned long int seed;
	seed = atoi(argv[2]);
	cout << "setting seed to " << seed << endl;
	gsl_rng_set(BGA::rng, seed);
    }
    if (argc > 3 && strlen(argv[3]) > 0)
    {
	snapshot_file = argv[3];
	cout << "molecule snapshots go to " << snapshot_file << endl;
    }

    // set lastMBadness to a maximum double
    numeric_limits<double> double_info;
    valarray<double> lastMBadness(double_info.max(), dtab->NAtoms+1);
    double best_largest = double_info.max();
    valarray<int> improved(1, logsize);

    Molecule mol(*dtab);
    mol.tol_dd = tol_dd;
    int fileno = 0;

    int maxatoms = 0;
    bool go_all_way = false;
    for (int trial = 0; ; ++trial)
    {
	// calculate pe
	double pe, impr_rate;
	if (mol.NAtoms() == mol.max_NAtoms())
	{
	    pe = 0.0;
	    go_all_way = false;
	}
	else if (mol.NAtoms() == 0)
	    pe = 1.0;
	else if (go_all_way)
	{
	    pe = 1.0;
	    if (mol.Badness() > 2*mol.NAtoms())
		go_all_way = false;
	}
	else
	{
	    impr_rate = 1.0*improved.sum()/improved.size();
	    pe = impr_rate*(pemax-pemin)+pemin;
	    if (impr_rate >= 0.66 && pallway > gsl_rng_uniform(BGA::rng))
		go_all_way = true;
	}
	cout << trial;
	if (pe > gsl_rng_uniform(BGA::rng))
	{
	    mol.Evolve(10,20,1000);
	    cout << "  Evolve()" << "  NAtoms = " << mol.NAtoms() << endl;
	}
	else
	{
	    int Npop = 0;
	    Npop = 1 + (int) floor(mol.Badness());
	    Npop = min(Npop, 5);
	    if (Npop > 1)
		Npop = 1 + gsl_rng_uniform_int(BGA::rng, Npop-1);
	    mol.Degenerate(Npop);
	    cout << "  Degenerate(" << Npop << ")  NAtoms = " << mol.NAtoms() <<  endl;
	}
	mol.PrintBadness();
	if (mol.NAtoms() == mol.max_NAtoms())
	    cout << "mol.Badness()/NAtoms = " << mol.Badness()/mol.max_NAtoms() << endl;
	// update lastMBadness and improved
	int ilog = trial % logsize;
	if (mol.Badness() < lastMBadness[mol.NAtoms()])
	{
	    if (mol.NAtoms() > maxatoms)
	    {
		best_largest = lastMBadness[mol.NAtoms()];
		maxatoms = mol.NAtoms();
	    }
	    lastMBadness[mol.NAtoms()] = mol.Badness();
	    improved[ilog] = 1;
	    /*
	    if (snapshot_file != NULL && mol.NAtoms() == maxatoms &&
		    best_largest >= mol.Badness())
	    {
		best_largest = mol.Badness();
		cout << "saving best molecule" << endl;
		mol.WriteXYZ(snapshot_file);
	    }
	    */
	    if (snapshot_file != NULL && improved[ilog])
	    {
		cout << "saving best molecule" << endl;
		char fname[255];
		sprintf(fname, "%s%04i", snapshot_file, fileno);
		sprintf(fname, "%s", snapshot_file);
		mol.WriteAtomEye(fname);
	    }
	    if (mol.NAtoms() == mol.max_NAtoms())
	    {
		if (mol.Badness() < avgmb*mol.max_NAtoms())
		{
		    cout << "that is solution!" << endl;
		    break;
		}
	    }
	}
	else
	{
	    improved[ilog] = 0;
	    if (lastMBadness[mol.NAtoms()] < 1e-4)
		lastMBadness[mol.NAtoms()] = 1e-4;
	}
	cout << endl;
    }
    return EXIT_SUCCESS;
}
