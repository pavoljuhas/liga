/***********************************************************************
* Short Title: MC molecule reconstruction from distance table
*
* Comments: Metropolis algorithm
*
* $Id$
***********************************************************************/

#include <limits>
#include <unistd.h>
#include <signal.h>
#include "ParseArgs.hpp"
#include "BGAlib.hpp"

////////////////////////////////////////////////////////////////////////
// RunPar_t
////////////////////////////////////////////////////////////////////////

struct RunPar_t
{
    RunPar_t();
    Molecule ProcessArguments(int argc, char * argv[]);
    // IO parameters
    string distfile;
    string inistru;
    string outstru;
    int saverate;
    int lograte;
    // MC parameters
    double boxsize;
    double kbt;
    double tol_bad;
    int seed;
    string penalty;
private:
    void print_help(ParseArgs& a);
    string version_string(string quote = "");
    list<string> validpars;
};

RunPar_t::RunPar_t()
{
    char *pnames[] = {
	"distfile", "inistru",
	"outstru", "saverate", "lograte", 
	"boxsize", "kbt", "tol_bad", "seed",
	"penalty" };
    validpars.insert(validpars.end(),
	    pnames, pnames+sizeof(pnames)/sizeof(char*));
}

void RunPar_t::print_help(ParseArgs& a)
{
    // /usage:/;/;/-s/.*/"&\\n"/
    // /cou/;/;/s/^\s*"\(.*\)\\n"/\1/ | '[put! ='/*' | /;/put ='*/'
    cout << 
"usage: " << a.cmd_t << "[-p PAR_FILE] [DISTFILE] [INISTRU] [par1=val1...]\n"
"run MC molecule simulation using distances from DISTFILE.  Parameters can\n"
"be set in PAR_FILE or on the command line, which overrides PAR_FILE.\n"
"Options:\n"
"  -p, --parfile=FILE    read parameters from FILE\n"
"  -h, --help            display this message\n"
"  -v, --version         show program version\n"
"IO parameters:\n"
"  distfile=FILE         target distance table - required\n"
"  inistru=FILE          initial structure - required\n"
"  outstru=FILE          where to save the best full molecule\n"
"  saverate=int          [10] minimum steps between outstru updates\n"
"  lograte=int           [100] minimum steps between log printout\n"
"Liga parameters\n"
"  boxsize=double        [0.1] size of box of possible MC step\n"
"  kbt=double            [0.001] Boltzman factor\n"
"  tol_bad=double        [1E-4] target normalized molecule badness\n"
"  seed=int              seed of random number generator\n"
"  penalty=string        dd penalty function [pow2], fabs, well\n"
;
}

string RunPar_t::version_string(string quote)
{
    using namespace std;
    ostringstream oss;
    oss << quote
        << "$Id$" << endl
#   if defined(__DATE__) && defined(__TIME__)
	<< quote << "compiled " __DATE__ " " __TIME__ << endl
#   endif
        ;
    return oss.str();
}

Molecule RunPar_t::ProcessArguments(int argc, char *argv[])
{
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
        exit(EXIT_FAILURE);
    }
    if (a.isopt("h") || argc == 1)
    {
        print_help(a);
        exit(EXIT_SUCCESS);
    }
    else if (a.isopt("v"))
    {
	cout << version_string();
        exit(EXIT_SUCCESS);
    }
    if (a.isopt("p"))
    {
        try {
            a.ReadPars(a.opts["p"].c_str());
        }
        catch (IOError(e)) {
            cerr << e.what() << endl;
            exit(EXIT_FAILURE);
        }
        catch (ParseArgsError(e)) {
            cerr << "invalid syntax in parameter file" << endl;
            cerr << e.what() << endl;
            exit(EXIT_FAILURE);
        }
    }
    try {
	a.ValidatePars(validpars);
    }
    catch (ParseArgsError(e)) {
	cerr << e.what() << endl;
	exit(EXIT_FAILURE);
    }
    // assign required parameters, distfile and inistru
    if (a.args.size() > 0)
        a.pars["distfile"] = a.args[0];
    if (a.args.size() > 1)
        a.pars["inistru"] = a.args[1];
    if (a.args.size() > 2)
    {
	cerr << argv[0] << ": too many file arguments" << endl;
	exit(EXIT_FAILURE);
    }
    if (!a.ispar("distfile"))
    {
        cerr << "Distance file not defined" << endl;
        exit(EXIT_FAILURE);
    }
    if (!a.ispar("inistru"))
    {
        cerr << "Initial structure not defined" << endl;
        exit(EXIT_FAILURE);
    }
    // intro messages
    string hashsep(72, '#');
    cout << hashsep << endl;
    cout << "# " << a.cmd_t << endl;
    cout << version_string("# ");
    char hostname[255];
    gethostname(hostname, 255);
    cout << "# " << hostname << endl;
    time_t cur_time = time(NULL);
    cout << "# " << ctime(&cur_time);
    cout << hashsep << endl;
    // distfile
    distfile = a.pars["distfile"];
    cout << "distfile=" << distfile << endl;
    DistanceTable* dtab;
    try {
        dtab = new DistanceTable(distfile.c_str());
    }
    catch (IOError) {
        exit(EXIT_FAILURE);
    }
    Molecule mol(*dtab);
    // inistru
    inistru = a.pars["inistru"];
    cout << "inistru=" << inistru << endl;
    try {
	mol.ReadXYZ(inistru.c_str());
    }
    catch (IOError) {
	exit(EXIT_FAILURE);
    }
    // outstru
    if (a.ispar("outstru"))
    {
        outstru = a.pars["outstru"];
        cout << "outstru=" << outstru << endl;
        saverate = a.GetPar<int>("saverate", 10);
        cout << "saverate=" << saverate << endl;
    }
    //  lograte
    lograte = a.GetPar<int>("lograte", 100);
    cout << "lograte=" << lograte << endl;
    // MC parameters
    try {
	// boxsize
	boxsize = a.GetPar<double>("boxsize", 0.001);
	cout << "boxsize=" << boxsize << endl;
	// kbt
	kbt = a.GetPar<double>("kbt", 0.1);
	cout << "kbt=" << kbt << endl;
	// tol_bad
	tol_bad = a.GetPar<double>("tol_bad", 1.0e-4);
	cout << "tol_bad=" << tol_bad << endl;
	mol.tol_nbad = tol_bad;
	// seed
	seed = a.GetPar<int>("seed", 0);
	if (seed)
	{
	    gsl_rng_set(BGA::rng, seed);
	    cout << "seed=" << seed << endl;
	}
	// penalty
	penalty = a.GetPar<string>("penalty", "pow2");
	if (penalty == "pow2")
	    mol.penalty = BGA::pow2;
	else if (penalty == "well")
	    mol.penalty = BGA::well;
	else if (penalty == "fabs")
	    mol.penalty = fabs;
	else
	{
	    cerr << "Invalid value of penalty parameter" << endl;
	    exit(EXIT_FAILURE);
	}
	cout << "penalty=" << penalty << endl;
    }
    catch (ParseArgsError(e)) {
	cerr << e.what() << endl;
	exit(EXIT_FAILURE);
    }
    // done
    cout << hashsep << endl << endl;
    return mol;
}


////////////////////////////////////////////////////////////////////////
// RunVar_t
////////////////////////////////////////////////////////////////////////

struct RunVar_t
{
    RunVar_t() : totsteps(0), accsteps(0), exiting(false)
    { }
    int totsteps;
    int accsteps;
    bool exiting;
};

////////////////////////////////////////////////////////////////////////
// SIGHUP handling
////////////////////////////////////////////////////////////////////////

int SIGHUP_received = 0;
void SIGHUP_handler(int signum)
{
    // die on 2nd call
    if (SIGHUP_received)
	exit(128 + signum);
    SIGHUP_received = signum;
}

////////////////////////////////////////////////////////////////////////
// Output helpers
////////////////////////////////////////////////////////////////////////

void save_outstru(Molecule& mol, RunPar_t& rp, RunVar_t& rv)
{
    static int cnt = 0;
    static double bestMNB = numeric_limits<double>().max();
    if (  rp.outstru.size() == 0 || rp.saverate == 0 ||
	    (++cnt < rp.saverate && !rv.exiting) )
	return;
    if (mol.NormBadness() < bestMNB)
    {
	bestMNB = mol.NormBadness();
	mol.WriteAtomEye(rp.outstru.c_str());
	cnt = 0;
    }
}

////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    // process arguments
    RunPar_t rp;
    RunVar_t rv;
    Molecule mol = rp.ProcessArguments(argc, argv);

    // print initial badness
    cout << rv.totsteps << ' ' << rv.accsteps << " I "
	<< mol.NormBadness() << endl;
    // store the best molecule ever
    Molecule best_mol(mol);
    cout << rv.totsteps << ' ' << rv.accsteps << " BC "
	<< best_mol.NormBadness() << endl;
    // watch for HUP
    signal(SIGHUP, SIGHUP_handler);
    // let the roulette begin
    while ( !(best_mol.NormBadness() < rp.tol_bad) )
    {
	if (SIGHUP_received)
	    break;
        ++rv.totsteps;
	// store original normalized badnees
	double nb0 = mol.NormBadness();
	// pick one atom
	int aidx = gsl_rng_uniform_int(BGA::rng, mol.NAtoms());
	Atom_t a0 = mol.Atom(aidx);
	// apply random shift to a copy
	double r1[3];
	r1[0] = a0.r[0] + rp.boxsize * (gsl_rng_uniform(BGA::rng) - 0.5);
	r1[1] = a0.r[1] + rp.boxsize * (gsl_rng_uniform(BGA::rng) - 0.5);
	r1[2] = a0.r[2] + rp.boxsize * (gsl_rng_uniform(BGA::rng) - 0.5);
	Atom_t a1(r1);
	// replace original atom with a new copy
	mol.Pop(aidx).Add(a1);
	double nb1 = mol.NormBadness();
	// accept according to Metropolis algorithm
	if (nb1 < nb0 || gsl_rng_uniform(BGA::rng) < exp(-(nb1-nb0)/rp.kbt) )
	    ++rv.accsteps;
	else
	{
	    // revert back, a1 is the last atom
	    mol.Pop(mol.NAtoms()-1).Add(a0);
	}
	// saving, log printing
        save_outstru(mol, rp, rv);
	if (rp.lograte && rv.totsteps % rp.lograte == 0)
	{
	    cout << rv.totsteps << ' ' << rv.accsteps << " L "
		<< mol.NormBadness() << endl;
	}
	if (nb1 < best_mol.NormBadness())
	{
	    best_mol = mol;
	    cout << rv.totsteps << ' ' << rv.accsteps << " BC "
		<< best_mol.NormBadness() << endl;
	}
    }
    cout << endl;
    if (SIGHUP_received)
	cout << "Received SIGHUP, graceful death." << endl << endl;
    else
	cout << "Solution found!!!" << endl << endl;
    BGA::cnt.PrintRunStats();
    // save last molecule
    rv.exiting = true;
    save_outstru(best_mol, rp, rv);
    if (SIGHUP_received)
	exit(SIGHUP+128);
    return EXIT_SUCCESS;
}
