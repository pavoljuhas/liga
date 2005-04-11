/***********************************************************************
* Short Title: molecule reconstruction from distance table
*
* Comments: self-tuning competition search for molecular configuration
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
    int natoms;
    string outstru;
    int saverate;
    string frames;
    int framesrate;
    // Walk parameters
    double tol_dd;
    double tol_bad;
    double maxcputime;
    int seed;
    double evolve_frac;
    bool evolve_relax;
    bool degenerate_relax;
    int ligasize;
    double stopgame;
    string penalty;
    int dist_trials;
    int tri_trials;
    int pyr_trials;
private:
    void print_help(ParseArgs& a);
    string version_string(string quote = "");
    list<string> validpars;
};

RunPar_t::RunPar_t()
{
    char *pnames[] = {
	"distfile", "inistru", "natoms",
	"outstru", "saverate", "frames", "framesrate",
	"tol_dd", "tol_bad", "maxcputime", "seed",
	"evolve_frac", "evolve_relax", "degenerate_relax",
	"ligasize", "stopgame",
	"penalty", "dist_trials", "tri_trials", "pyr_trials" };
    validpars.insert(validpars.end(),
	    pnames, pnames+sizeof(pnames)/sizeof(char*));
}

void RunPar_t::print_help(ParseArgs& a)
{
    // /usage:/;/;/-s/.*/"&\\n"/
    // /cou/;/;/s/^\s*"\(.*\)\\n"/\1/ | '[put! ='/*' | /;/put ='*/'
    cout << 
"usage: " << a.cmd_t << " [-p PAR_FILE] [DISTFILE] [par1=val1 par2=val2...]\n"
"run gizaliga simulation using distances from DISTFILE.  Parameters can\n"
"be set in PAR_FILE or on the command line, which overrides PAR_FILE.\n"
"Options:\n"
"  -p, --parfile=FILE    read parameters from FILE\n"
"  -h, --help            display this message\n"
"  -v, --version         show program version\n"
"IO parameters:\n"
"  distfile=FILE         target distance table\n"
"  inistru=FILE          initial structure [empty box]\n"
"  natoms=int            size of molecule, use with loose target distances\n"
"  outstru=FILE          where to save the best full molecule\n"
"  saverate=int          [10] minimum iterations between outstru updates\n"
"  frames=FILE           save intermediate structures to FILE.liga_round\n"
"  framesrate=int        [10] number of iterations between frame saves\n"
"Liga parameters\n"
"  tol_dd=double         [0.1] distance is not used when dd=|d-d0|>=tol_dd\n"
"  tol_bad=double        [1E-4] target normalized molecule badness\n"
"  maxcputime=double     [0] when set, maximum CPU time in seconds\n"
"  seed=int              seed of random number generator\n"
"  evolve_frac=double    [0.1] fraction of tol_bad threshold of tested atoms\n"
"  evolve_relax=bool     [false] relax the worst atom after addition\n"
"  degenerate_relax=bool [false] relax the worst atom after removal\n"
"  ligasize=int          [10] number of teams per division\n"
"  stopgame=double       [0.025] skip division when winner is worse\n"
"  penalty=string        dd penalty function [pow2], fabs, well\n"
"  dist_trials=int       [10] good distance atoms to try\n"
"  tri_trials=int        [20] godd triangle atoms to try\n"
"  pyr_trials=int        [1000] good pyramid atoms to try\n"
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
    // assign run parameters
    // distfile
    if (a.args.size())
        a.pars["distfile"] = a.args[0];
    if (a.args.size() > 1)
    {
	cerr << argv[0] << ": several DISTFILE arguments" << endl;
	exit(EXIT_FAILURE);
    }
    if (!a.ispar("distfile"))
    {
        cerr << "Distance file not defined" << endl;
        exit(EXIT_FAILURE);
    }
    distfile = a.pars["distfile"];
    DistanceTable* dtab;
    try {
        dtab = new DistanceTable(distfile.c_str());
    }
    catch (IOError) {
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
    Molecule mol(*dtab);
    cout << "distfile=" << distfile << endl;
    // inistru
    if (a.ispar("inistru"))
    {
	inistru = a.pars["inistru"];
	cout << "inistru=" << inistru << endl;
	try {
	    mol.ReadXYZ(inistru.c_str());
	}
	catch (IOError) {
	    exit(EXIT_FAILURE);
	}
    }
    // natoms
    if (a.ispar("natoms"))
    {
	try {
	    natoms = a.GetPar<int>("natoms");
	    cout << "natoms=" << natoms << endl;
	    mol.Set_max_NAtoms(natoms);
	}
        catch (ParseArgsError(e)) {
            cerr << e.what() << endl;
            exit(EXIT_FAILURE);
	}
	catch (InvalidMolecule) {
	    exit(EXIT_FAILURE);
	}
    }
    // outstru
    if (a.ispar("outstru"))
    {
        outstru = a.pars["outstru"];
        cout << "outstru=" << outstru << endl;
        saverate = a.GetPar<int>("saverate", 10);
        cout << "saverate=" << saverate << endl;
    }
    // frames, framesrate
    if (a.ispar("frames"))
    {
        frames = a.pars["frames"];
        cout << "frames=" << frames << endl;
        framesrate = a.GetPar<int>("framesrate", 10);
        cout << "framesrate=" << framesrate << endl;
    }
    // liga parameters
    try {
	// tol_dd
	tol_dd = a.GetPar<double>("tol_dd", 0.1);
	cout << "tol_dd=" << tol_dd << endl;
	mol.tol_dd = tol_dd;
	// tol_bad
	tol_bad = a.GetPar<double>("tol_bad", 1.0e-4);
	cout << "tol_bad=" << tol_bad << endl;
	mol.tol_nbad = tol_bad;
	// maxcputime
	maxcputime = a.GetPar<double>("maxcputime", 0.0);
	if (maxcputime > 0.0)
	    cout << "maxcputime=" << maxcputime << endl;
	// seed
	seed = a.GetPar<int>("seed", 0);
	if (seed)
	{
	    gsl_rng_set(BGA::rng, seed);
	    cout << "seed=" << seed << endl;
	}
	// evolve_frac
	evolve_frac = a.GetPar<double>("evolve_frac", 0.1);
	cout << "evolve_frac=" << evolve_frac << endl;
	mol.evolve_frac = evolve_frac;
	// evolve_relax
	evolve_relax = a.GetPar<bool>("evolve_relax", false);
	cout << "evolve_relax=" << evolve_relax << endl;
	mol.evolve_relax = evolve_relax;
	// degenerate_relax
	degenerate_relax = a.GetPar<bool>("degenerate_relax", false);
	cout << "degenerate_relax=" << degenerate_relax << endl;
	mol.degenerate_relax = degenerate_relax;
	// ligasize
	ligasize = a.GetPar<int>("ligasize", 10);
	cout << "ligasize=" << ligasize << endl;
	// stopgame
	stopgame = a.GetPar<double>("stopgame", 0.0025);
	cout << "stopgame=" << stopgame << endl;
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
	// dist_trials
	dist_trials = a.GetPar("dist_trials", 10);
	cout << "dist_trials=" << dist_trials << endl;
	// tri_trials
	tri_trials = a.GetPar("tri_trials", 20);
	cout << "tri_trials=" << tri_trials << endl;
	// pyr_trials
	pyr_trials = a.GetPar("pyr_trials", 1000);
	cout << "pyr_trials=" << pyr_trials << endl;
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
    RunVar_t() : liga_round(0), full_liga(false), exiting(false)
    { }
    int liga_round;
    bool full_liga;
    bool exiting;
};


////////////////////////////////////////////////////////////////////////
// Division_t
////////////////////////////////////////////////////////////////////////

typedef Molecule* PMOL;
struct Division_t : public vector<Molecule*>
{
public:
    // constructors
    Division_t(int s)  : vector<PMOL>(), max_size(s) { }
    Division_t(const Division_t& div0) :
	vector<PMOL>(div0), max_size(div0.max_size) { }
    ~Division_t()
    {
	for (iterator ii = begin(); ii != end(); ++ii)
	    delete *ii;
    }
    Division_t& operator= (const vector<PMOL>& div0)
    {
	*this = div0;
    }
    Division_t& operator= (const Division_t& div0)
    {
	*this = vector<PMOL>(div0);
	max_size = div0.max_size;
    }
    int find_winner();
    PMOL& winner();
    int find_looser();
    PMOL& looser();
    PMOL& best();
    inline bool full() { return !(size() < max_size); }
    inline int fullsize() { return max_size; }
private:
    mutable int max_size;
};

int Division_t::find_winner()
{
    // evaluate molecule fitness
    valarray<double> vmfit(size());
    double *pd = &vmfit[0];
    for (iterator mi = begin(); mi != end(); ++mi, ++pd)
	*pd = (*mi)->NormBadness();
    // then get the reciprocal value
    vmfit = vdrecipw0(vmfit);
    double *mfit = &vmfit[0];
    int idx = *(random_wt_choose(1, mfit, size()).begin());
    return idx;
}

PMOL& Division_t::winner()
{
    return at(find_winner());
}

int Division_t::find_looser()
{
    // evaluate molecule fitness
    valarray<double> vmbad(size());
    double *pd = &vmbad[0];
    for (iterator mi = begin(); mi != end(); ++mi, ++pd)
	*pd = (*mi)->NormBadness();
    double *mbad = &vmbad[0];
    int idx = *(random_wt_choose(1, mbad, size()).begin());
    return idx;
}

PMOL& Division_t::looser()
{
    return at(find_looser());
}

bool comp_PMOL_Badness(const PMOL& lhs, const PMOL& rhs)
{
    return lhs->Badness() < rhs->Badness();
}

PMOL& Division_t::best()
{
    // evaluate molecule fitness
    iterator pm = min_element(begin(), end(), comp_PMOL_Badness);
    return *pm;
}


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
    double dbmax = numeric_limits<double>().max();
    static valarray<double> bestMNB(dbmax, mol.max_NAtoms()+1);
    if (  rp.outstru.size() == 0 || rp.saverate == 0 ||
	    (++cnt < rp.saverate && !rv.exiting) )
	return;
    if (mol.NormBadness() < bestMNB[mol.NAtoms()])
    {
	bestMNB[mol.NAtoms()] = mol.NormBadness();
	mol.WriteAtomEye(rp.outstru.c_str());
	cnt = 0;
    }
}

void save_frames(Molecule& mol, RunPar_t& rp, RunVar_t& rv)
{
    static int cnt = 0;
    DistanceTable dummydtgt;
    static Molecule last_frame(dummydtgt);
    if (  rp.frames.size() == 0 || rp.framesrate == 0 ||
	    (++cnt < rp.framesrate && !rv.exiting) ||
	    last_frame == mol
       )
	return;
    ostringstream oss;
    oss << rp.frames << "." << rv.liga_round;
    mol.WriteAtomEye(oss.str().c_str());
    last_frame = mol;
    cnt = 0;
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

    // initialize liga divisions, primitive divisions have only 1 team
    vector<Division_t> liga;
    for (int i = 0; i <= mol.max_NAtoms(); ++i)
    {
        int divsize = (i < 2) ? 1 : rp.ligasize;
        liga.push_back(Division_t(divsize));
    }
    // put initial molecule to its division
    PMOL first_team = new Molecule(mol);
    cout << "Initial team" << endl;
    cout << rv.liga_round << " I " << first_team->NAtoms() << ' '
	<< first_team->NormBadness() << endl;
    liga[mol.NAtoms()].push_back(first_team);
    // fill lower divisions
    cout << "Filling lower divisions" << endl;
    for (int level = mol.NAtoms()-1; level >= 0; --level)
    {
        PMOL parent_team = liga[level+1].back();
        PMOL lower_team = new Molecule(*parent_team);
        lower_team->Degenerate(1);
	cout << rv.liga_round << " L " << lower_team->NAtoms() << ' '
	    << lower_team->NormBadness() << endl;
        liga[level].push_back(lower_team);
    }
    cout << "Done" << endl;
    // the first world champion is the initial molecule
    PMOL world_champ = first_team;
    cout << rv.liga_round << " WC " << world_champ->NAtoms() << ' '
	<< world_champ->NormBadness() << endl;
    // save the best team ever:
    Molecule best_champ(*world_champ);
    cout << rv.liga_round << " BC " << best_champ.NAtoms() << ' '
	<< best_champ.NormBadness() << endl;
    cout << endl << "Starting the game ... now." << endl;
    // watch for HUP
    signal(SIGHUP, SIGHUP_handler);
    // let the game begin
    while ( !(world_champ->Full() && world_champ->NormBadness() < rp.tol_bad) )
    {
	// check special stoping conditions
	if (SIGHUP_received)
	    break;
	if (rp.maxcputime > 0.0 && BGA::cnt.CPUTime() > rp.maxcputime)
	    break;
        ++rv.liga_round;
	typedef vector<Division_t>::iterator VDit;
	int lo_level = 0;
	for (VDit lo_div = liga.begin();
		lo_div < liga.end()-1; ++lo_div, ++lo_level)
	{
	    if (lo_div->size() == 0)
		continue;
	    int winner_idx = lo_div->find_winner();
	    PMOL advancing = lo_div->at(winner_idx);
	    if (! (advancing->NormBadness() < rp.stopgame) )
		continue;
	    double adv_bad0 = advancing->NormBadness();
	    if (! lo_div->full() )
	    {
		// save clone of advancing winner
		PMOL winner_clone = new Molecule(*advancing);
		lo_div->push_back(winner_clone);
	    }
	    // advance as far as possible
	    advancing->Evolve(rp.dist_trials, rp.tri_trials, rp.pyr_trials);
	    int hi_level = advancing->NAtoms();
	    VDit hi_div = liga.begin() + hi_level;
	    if (hi_div->size() == 0)
		hi_div->push_back(new Molecule(*advancing));
	    int looser_idx = hi_div->find_looser();
	    PMOL descending = hi_div->at(looser_idx);
	    double desc_bad0 = descending->NormBadness();
	    if (! hi_div->full() )
	    {
		// save clone of descending looser
		PMOL looser_clone = new Molecule(*descending);
		hi_div->push_back(looser_clone);
	    }
	    // clone winner if he made a good advance
	    if (descending->NormBadness() > advancing->NormBadness())
	    {
		*descending = *advancing;
	    }
	    descending->Degenerate(hi_level-lo_level);
	    // all set now so we can swap winner and looser
	    (*hi_div)[looser_idx] = advancing;
	    (*lo_div)[winner_idx] = descending;
	    cout << rv.liga_round;
	    cout << " A " <<
		lo_level << ' ' << adv_bad0 << ' ' <<
		hi_level << ' ' << advancing->NormBadness() << "    ";
	    cout << " D " <<
		hi_level << ' ' << desc_bad0 << ' ' <<
		lo_level << ' ' << descending->NormBadness() << endl;
	    // no need to finish round if we found champion
	    if (advancing->Full() && advancing->NormBadness() < rp.tol_bad)
		break;
	    if (SIGHUP_received)
		break;
	}
	if (liga.back().size() != 0)
	    world_champ = liga.back().best();
	else
	{
	    typedef vector<Division_t>::reverse_iterator VDrit;
	    for (VDrit ii = liga.rbegin(); ii != liga.rend(); ++ii)
	    {
		if (ii->size())
		{
		    world_champ = ii->best();
		    break;
		}
	    }
	}
	cout << rv.liga_round << " WC " << world_champ->NAtoms() << ' '
	    << world_champ->NormBadness() << endl;
	if (    world_champ->NAtoms() > best_champ.NAtoms() ||
		world_champ->NormBadness() < best_champ.NormBadness()
	   )
	{
	    best_champ = *world_champ;
	    cout << rv.liga_round << " BC " << best_champ.NAtoms() << ' '
		<< best_champ.NormBadness() << endl;
	}
        save_outstru(best_champ, rp, rv);
        save_frames(best_champ, rp, rv);
    }
    cout << endl;
    int exit_code;
    if (SIGHUP_received)
    {
	cout << "Received SIGHUP, graceful death." << endl << endl;
	exit_code = SIGHUP+128;
    }
    else if (rp.maxcputime > 0.0 && BGA::cnt.CPUTime() > rp.maxcputime)
    {
	cout << "Exceeded maxcputime." << endl << endl;
	exit_code = 1;
    }
    else
    {
	cout << "Solution found!!!" << endl << endl;
	exit_code = EXIT_SUCCESS;
    }
    BGA::cnt.PrintRunStats();
    // save last frame
    rv.exiting = true;
    save_outstru(best_champ, rp, rv);
    save_frames(best_champ, rp, rv);
    if (SIGHUP_received)
	exit(exit_code);
    return exit_code;
}
