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
// helper types: TraceId_t
////////////////////////////////////////////////////////////////////////

struct TraceId_t
{
    int season;
    int ini_level;
    int fin_level;
};

const int EXIT_INPUT_ERROR = 2;

////////////////////////////////////////////////////////////////////////
// RunPar_t
////////////////////////////////////////////////////////////////////////

struct RunPar_t
{
    RunPar_t();
    Molecule ProcessArguments(int argc, char * argv[]);
    // Output option
    bool quiet;
    bool trace;
    // IO parameters
    string distfile;
    string inistru;
    string outstru;
    string outfmt;
    int saverate;
    bool saveall;
    string frames;
    int framesrate;
    vector<TraceId_t> framestrace;
    // Liga parameters
    double tol_dd;
    double tol_bad;
    int natoms;
    vector<int> fixed_atoms;
    int centersize;
    double maxcputime;
    int seed;
    double evolve_frac;
    bool evolve_relax;
    bool degenerate_relax;
    int ligasize;
    double stopgame;
    int dist_trials;
    int tri_trials;
    int pyr_trials;
    double lookout_prob;
    // Constrains
    vector<double> bangle_range;
    double max_dist;

private:
    void print_help(ParseArgs& a);
    string version_string(string quote = "");
    void print_pars(ParseArgs& a);
    list<string> validpars;
};

RunPar_t::RunPar_t()
{
    char *pnames[] = {
	"distfile", "inistru",
	"outstru", "outfmt", "saverate", "saveall",
	"frames", "framesrate", "framestrace",
	"tol_dd", "tol_bad", "natoms", "fixed_atoms", "centersize",
	"maxcputime", "seed",
	"evolve_frac", "evolve_relax", "degenerate_relax",
	"ligasize", "stopgame",
	"dist_trials", "tri_trials", "pyr_trials", "lookout_prob",
	"bangle_range", "max_dist" };
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
"  -q, --quiet           restrict output to champion quality\n"
"  -t, --trace           keep and show frame trace of the champion\n"
"  -h, --help            display this message\n"
"  -V, --version         show program version\n"
"IO parameters:\n"
"  distfile=FILE         target distance table\n"
"  inistru=FILE          initial structure [empty box]\n"
"  outstru=FILE          where to save the best full molecule\n"
"  outfmt=string         [xyz], atomeye - outstru file format\n"
"  saverate=int          [10] minimum iterations between outstru updates\n"
"  saveall=bool          [false] save best molecules from all divisions\n"
"  frames=FILE           save intermediate structures to FILE.season\n"
"  framesrate=int        [0] number of iterations between frame saves\n"
"  framestrace=array     [] triples of (season, initial, final level)\n"
"Liga parameters:\n"
"  tol_dd=double         [0.1] distance is not used when dd=|d-d0|>=tol_dd\n"
"  tol_bad=double        [1E-4] target normalized molecule badness\n"
"  natoms=int            use with loose distfiles or when tol_dd==0\n"
"  fixed_atoms=ranges    [] indices of fixed atoms in inistru (start at 1)\n"
"  centersize=int        [0] shift smaller molecules to the origin\n"
"  maxcputime=double     [0] when set, maximum CPU time in seconds\n"
"  seed=int              seed of random number generator\n"
"  evolve_frac=double    [0.1] fraction of tol_bad threshold of tested atoms\n"
"  evolve_relax=bool     [false] relax the worst atom after addition\n"
"  degenerate_relax=bool [false] relax the worst atom after removal\n"
"  ligasize=int          [10] number of teams per division\n"
"  stopgame=double       [0.0025] skip division when winner is worse\n"
"  dist_trials=int       [10] good distance atoms to try\n"
"  tri_trials=int        [20] godd triangle atoms to try\n"
"  pyr_trials=int        [1000] good pyramid atoms to try\n"
"  lookout_prob=double   [0.0] lookout probability for 2nd and 3rd atoms\n"
"Constrains (applied only when set):\n"
"  bangle_range=array    (max_blen, low[, high]) bond angle constraint\n"
"  max_dist=double       distance limit for rejecting lone atoms\n"
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

void RunPar_t::print_pars(ParseArgs& a)
{
    // print out all run parameters
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
    cout << "distfile=" << distfile << endl;
    // inistru
    if (a.ispar("inistru"))
    {
	cout << "inistru=" << inistru << endl;
    }
    // outstru, outfmt, saverate, saveall
    if (a.ispar("outstru"))
    {
        cout << "outstru=" << outstru << endl;
	cout << "outfmt=" << outfmt << endl;
        cout << "saverate=" << saverate << endl;
	cout << "saveall=" << saveall << endl;
    }
    // frames, framestrace, framesrate
    if (a.ispar("frames"))
    {
	cout << "frames=" << frames << endl;
	// framestrace
	if (framestrace.size() != 0)
	{
	    cout << "framestrace=";
	    for (vector<TraceId_t>::iterator ii = framestrace.begin();
		    ii != framestrace.end(); ++ii)
	    {
		cout << " \\" << endl;
		cout << "    " << ii->season;
		cout << ' ' << ii->ini_level;
		cout << ' ' << ii->fin_level;
	    }
	    cout << endl;
	}
	// framesrate - print only if framestrace is empty
	else
	{
	    cout << "framesrate=" << framesrate << endl;
	}
    }
    // liga parameters
    // tol_dd, tol_bad 
    cout << "tol_dd=" << tol_dd << endl;
    cout << "tol_bad=" << tol_bad << endl;
    // natoms
    if (a.ispar("natoms"))
    {
	cout << "natoms=" << natoms << endl;
    }
    // fixed_atoms
    if (a.ispar("fixed_atoms") && !fixed_atoms.empty())
    {
	cout << "fixed_atoms=";
	for (int start = 0; start < fixed_atoms.size(); )
	{
	    int range = 0;
	    for ( ; start + range < fixed_atoms.size() &&
		    fixed_atoms[start]+range == fixed_atoms[start+range];
		    ++range )
	    { }
	    if (start > 0)  cout << ',';
	    cout << fixed_atoms[start];
	    if (range > 1)
	    {
		cout << ".." << fixed_atoms[start+range-1];
	    }
	    start += range;
	}
	cout << endl;
    }
    // centersize
    if (a.ispar("centersize"))
    {
	cout << "centersize=" << centersize << endl;
    }
    // maxcputime
    if (maxcputime > 0.0)
    {
	cout << "maxcputime=" << maxcputime << endl;
    }
    // seed
    if (seed)
    {
	cout << "seed=" << seed << endl;
    }
    // evolve_frac, evolve_relax, degenerate_relax
    cout << "evolve_frac=" << evolve_frac << endl;
    cout << "evolve_relax=" << evolve_relax << endl;
    cout << "degenerate_relax=" << degenerate_relax << endl;
    // ligasize, stopgame, dist_trials, tri_trials, pyr_trials, lookout_prob
    cout << "ligasize=" << ligasize << endl;
    cout << "stopgame=" << stopgame << endl;
    cout << "dist_trials=" << dist_trials << endl;
    cout << "tri_trials=" << tri_trials << endl;
    cout << "pyr_trials=" << pyr_trials << endl;
    cout << "lookout_prob=" << lookout_prob << endl;
    // constraints
    // bangle_range
    if (a.ispar("bangle_range"))
    {
	cout << "bangle_range=" << bangle_range[0];
	for (int i = 1; i < bangle_range.size(); ++i)
	{
	    cout << ',' << bangle_range[i];
	}
	cout << endl;
    }
    // max_dist
    if (a.ispar("max_dist"))
    {
	cout << "max_dist=" << max_dist << endl;
    }
    // finish done
    cout << hashsep << endl << endl;
}

Molecule RunPar_t::ProcessArguments(int argc, char *argv[])
{
    char *short_options =
        "p:qthV";
    // parameters and options
    option long_options[] = {
        {"parfile", 1, 0, 'p'},
        {"quiet", 0, 0, 'q'},
        {"trace", 0, 0, 't'},
        {"help", 0, 0, 'h'},
        {"version", 0, 0, 'V'},
        {0, 0, 0, 0}
    };
    ParseArgs a(argc, argv, short_options, long_options);
    try {
        a.Parse();
    }
    catch (ParseArgsError) {
        exit(EXIT_INPUT_ERROR);
    }
    if (a.isopt("h") || argc == 1)
    {
        print_help(a);
        exit(EXIT_SUCCESS);
    }
    else if (a.isopt("V"))
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
            exit(EXIT_INPUT_ERROR);
        }
        catch (ParseArgsError(e)) {
            cerr << "invalid syntax in parameter file" << endl;
            cerr << e.what() << endl;
            exit(EXIT_INPUT_ERROR);
        }
    }
    quiet = a.isopt("q");
    trace = a.isopt("t");
    try {
	a.ValidatePars(validpars);
    }
    catch (ParseArgsError(e)) {
	cerr << e.what() << endl;
	exit(EXIT_INPUT_ERROR);
    }
    // assign run parameters
    // distfile
    if (a.args.size())
        a.pars["distfile"] = a.args[0];
    if (a.args.size() > 1)
    {
	cerr << argv[0] << ": several DISTFILE arguments" << endl;
	exit(EXIT_INPUT_ERROR);
    }
    if (!a.ispar("distfile"))
    {
        cerr << "Distance file not defined" << endl;
        exit(EXIT_INPUT_ERROR);
    }
    distfile = a.pars["distfile"];
    DistanceTable* dtab;
    try {
        dtab = new DistanceTable(distfile.c_str());
    }
    catch (IOError) {
        exit(EXIT_INPUT_ERROR);
    }
    // create empty molecule
    Molecule mol(*dtab);
    // inistru
    if (a.ispar("inistru"))
    {
	inistru = a.pars["inistru"];
	try {
	    mol.ReadXYZ(inistru.c_str());
	}
	catch (IOError) {
	    exit(EXIT_INPUT_ERROR);
	}
    }
    // outstru, outfmt, saverate, saveall
    if (a.ispar("outstru"))
    {
	// outstru
        outstru = a.pars["outstru"];
	// outfmt
	outfmt = a.GetPar<string>("outfmt", "xyz");
	if (outfmt == "xyz")
	    mol.OutFmtXYZ();
	else if (outfmt == "atomeye")
	    mol.OutFmtAtomEye();
	else
	{
	    cerr << "Invalid value of outfmt parameter" << endl;
	    exit(EXIT_INPUT_ERROR);
	}
	// saverate
        saverate = a.GetPar<int>("saverate", 10);
	// saveall
	saveall = a.GetPar<bool>("saveall", false);
    }
    // frames, framestrace, framesrate
    if (a.ispar("frames"))
    {
	// frames
	frames = a.pars["frames"];
	// framestrace
	if (a.ispar("framestrace"))
	{
	    framesrate = 0;
	    vector<int> stp = a.GetParVec<int>("framestrace");
	    if (stp.size() % 3)
	    {
		cerr << "framestrace must have 3n entries" << endl;
		exit(EXIT_INPUT_ERROR);
	    }
	    // check if seasons in stp are ordered
	    for (int i = 3; i < stp.size(); i += 3)
	    {
		if (stp[i] < stp[i-3])
		{
		    cerr << "framestrace seasons must be sorted" << endl;
		    exit(EXIT_INPUT_ERROR);
		}
	    }
	    TraceId_t tid;
	    for (int i = 0; i < stp.size(); i += 3)
	    {
		tid.season = stp[i];
		tid.ini_level = stp[i+1];
		tid.fin_level = stp[i+2];
		framestrace.push_back(tid);
	    }
	}
	// framesrate - parse only if framestrace is empty
	if (framestrace.size() == 0)
	{
	    framesrate = a.GetPar<int>("framesrate", 0);
	}
    }
    // liga parameters
    try {
	// tol_dd
	tol_dd = a.GetPar<double>("tol_dd", 0.1);
	mol.tol_dd = tol_dd;
	// tol_bad
	tol_bad = a.GetPar<double>("tol_bad", 1.0e-4);
	mol.tol_nbad = tol_bad;
	// natoms must be set after tol_dd
	if (a.ispar("natoms"))
	{
	    try {
		natoms = a.GetPar<int>("natoms");
		mol.Set_max_NAtoms(natoms);
	    }
	    catch (InvalidMolecule) {
		exit(EXIT_INPUT_ERROR);
	    }
	}
	// fixed_atoms must be set after inistru
	if (a.ispar("fixed_atoms"))
	{
	    fixed_atoms = a.ExpandRangePar("fixed_atoms");
	    sort(fixed_atoms.begin(), fixed_atoms.end());
	    vector<int>::iterator last;
	    last = unique(fixed_atoms.begin(), fixed_atoms.end());
	    fixed_atoms.erase(last, fixed_atoms.end());
	    vector<int>::const_iterator ii;
	    try {
		for (ii = fixed_atoms.begin(); ii != fixed_atoms.end(); ++ii)
		{
		    mol.Fix(*ii - 1);
		}
	    }
	    catch (range_error) {
		cerr << "fixed_atoms - invalid index " << *ii << endl;
		exit(EXIT_INPUT_ERROR);
	    }
	}
	// centersize
	centersize = a.GetPar<int>("centersize", 0);
	mol.center_size = centersize;
	// maxcputime
	maxcputime = a.GetPar<double>("maxcputime", 0.0);
	// seed
	seed = a.GetPar<int>("seed", 0);
	if (seed)
	{
	    gsl_rng_set(BGA::rng, seed);
	}
	// evolve_frac
	evolve_frac = a.GetPar<double>("evolve_frac", 0.1);
	mol.evolve_frac = evolve_frac;
	// evolve_relax
	evolve_relax = a.GetPar<bool>("evolve_relax", false);
	mol.evolve_relax = evolve_relax;
	// degenerate_relax
	degenerate_relax = a.GetPar<bool>("degenerate_relax", false);
	mol.degenerate_relax = degenerate_relax;
	// ligasize
	ligasize = a.GetPar<int>("ligasize", 10);
	// stopgame
	stopgame = a.GetPar<double>("stopgame", 0.0025);
	// dist_trials
	dist_trials = a.GetPar("dist_trials", 10);
	// tri_trials
	tri_trials = a.GetPar("tri_trials", 20);
	// pyr_trials
	pyr_trials = a.GetPar("pyr_trials", 1000);
	// lookout_prob
	lookout_prob = a.GetPar("lookout_prob", 0.0);
	// bangle_range
	if (a.ispar("bangle_range"))
	{
	    bangle_range = a.GetParVec<double>("bangle_range");
	    if (bangle_range.size() != 2 && bangle_range.size() != 3)
	    {
		cerr << "bangle_range must have 2 or 3 arguments" << endl;
		exit(EXIT_INPUT_ERROR);
	    }
	    double max_blen = bangle_range[0];
	    BondAngleFilter_t* pbaf = new BondAngleFilter_t(max_blen);
	    pbaf->lo_bangle = bangle_range[1];
	    if (bangle_range.size() == 3)
	    {
		pbaf->hi_bangle = bangle_range[2];
	    }
	    mol.atom_filters.push_back(pbaf);
	}
	// max_dist
	if (a.ispar("max_dist"))
	{
	    max_dist = a.GetPar<double>("max_dist");
	    LoneAtomFilter_t* plaf = new LoneAtomFilter_t(max_dist);
	    plaf->max_dist = max_dist;
	    mol.atom_filters.push_back(plaf);
	}
    }
    catch (ParseArgsError(e)) {
	cerr << e.what() << endl;
	exit(EXIT_INPUT_ERROR);
    }
    // done
    print_pars(a);
    return mol;
}


////////////////////////////////////////////////////////////////////////
// RunVar_t
////////////////////////////////////////////////////////////////////////

struct RunVar_t
{
    RunVar_t() : season(0), exiting(false)
    { }
    int season;
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
	return *this;
    }
    Division_t& operator= (const Division_t& div0)
    {
	*this = vector<PMOL>(div0);
	max_size = div0.max_size;
	return *this;
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
    int idx = random_wt_choose(1, mfit, size()).front();
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
    int idx = random_wt_choose(1, mbad, size()).front();
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

void save_outstru(vector<Division_t>& liga, RunPar_t& rp, RunVar_t& rv)
{
    static int cnt = 0;
    double dbmax = numeric_limits<double>().max();
    static valarray<double> bestMNB(dbmax, liga.size());
    if (  rp.outstru.size() == 0 || rp.saverate == 0 ||
	    (++cnt < rp.saverate && !rv.exiting) )
	return;
    // start saving from the largest non-empty division
    for (int level = liga.size() - 1; level > 0; --level)
    {
	if ( liga[level].empty() )
	    continue;
	PMOL level_champ = liga[level].best();
	// save only when there is clear improvement
	if ( eps_lt(level_champ->NormBadness(), bestMNB[level]) )
	    bestMNB[level] = level_champ->NormBadness();
	else if (rp.saveall)
	    continue;
	else
	    break;
	// here we have something to save
	cnt = 0;
	// obtain file name
	string fname = rp.outstru;
	if ( rp.saveall )
	{
	    ostringstream oss;
	    oss << ".L" << level_champ->NAtoms();
	    fname.append( oss.str() );
	}
	level_champ->WriteFile(fname.c_str());
	// jump out if we are not saving all divisions
	if ( ! rp.saveall )
	    break;
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
    oss << rp.frames << "." << rv.season;
    mol.WriteAtomEye(oss.str().c_str());
    last_frame = mol;
    cnt = 0;
}

// convert to Liga_t member
void save_frames_trace(PMOL advancing, PMOL descending, RunPar_t& rp,
	RunVar_t& rv)
{
    static vector<TraceId_t>::iterator ti = rp.framestrace.begin();
    if (ti == rp.framestrace.end())
	return;
    if (rv.season != ti->season)
	return;
    if (ti->ini_level < ti->fin_level && ti->fin_level == advancing->NAtoms())
    {
	int trace_no = int(ti - rp.framestrace.begin()) + 1;
	ostringstream oss;
	oss << rp.frames << "." << rv.season << '.' << trace_no;
	advancing->WriteAtomEye(oss.str().c_str());
	if (++ti == rp.framestrace.end())
	    return;
    }
    if (ti->ini_level > ti->fin_level && ti->fin_level == descending->NAtoms())
    {
	int trace_no = int(ti - rp.framestrace.begin()) + 1;
	ostringstream oss;
	oss << rp.frames << "." << rv.season << '.' << trace_no;
	descending->WriteAtomEye(oss.str().c_str());
	if (++ti == rp.framestrace.end())
	    return;
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
    int base_level = mol.NFixed();
    // initialize liga divisions, primitive divisions have only 1 team
    vector<Division_t> liga;
    for (int i = 0; i <= mol.max_NAtoms(); ++i)
    {
        int divsize = rp.ligasize;
	if (i < base_level)  divsize = 0;
	else if (i < 2)      divsize = 1;
        liga.push_back(Division_t(divsize));
    }
    // put initial molecule to its division
    PMOL first_team = new Molecule(mol);
    cout << "Initial team" << endl;
    cout << rv.season << " I " << first_team->NAtoms() << ' '
	<< first_team->NormBadness() << endl;
    liga[mol.NAtoms()].push_back(first_team);
    // fill lower divisions
    cout << "Filling lower divisions" << endl;
    for (int level = mol.NAtoms()-1; level >= base_level; --level)
    {
        PMOL parent_team = liga[level+1].back();
        PMOL lower_team = new Molecule(*parent_team);
        lower_team->Degenerate(1);
	cout << rv.season << " L " << lower_team->NAtoms() << ' '
	    << lower_team->NormBadness() << endl;
        liga[level].push_back(lower_team);
    }
    cout << "Done" << endl;
    // the first world champion is the initial molecule
    PMOL world_champ = first_team;
    cout << rv.season << " WC " << world_champ->NAtoms() << ' '
	<< world_champ->NormBadness() << endl;
    // save the best team ever:
    Molecule best_champ(*world_champ);
    cout << rv.season << " BC " << best_champ.NAtoms() << ' '
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
        ++rv.season;
	typedef vector<Division_t>::iterator VDit;
	VDit lo_div = liga.begin() + base_level;
	int lo_level = base_level;
	for ( ; lo_div < liga.end()-1; ++lo_div, ++lo_level)
	{
	    if (lo_div->size() == 0)
		continue;
	    int winner_idx = lo_div->find_winner();
	    PMOL advancing = lo_div->at(winner_idx);
	    if (! (advancing->NormBadness() < rp.stopgame) )
		continue;
	    bool advancing_best = (advancing == lo_div->best());
	    double adv_bad0 = advancing->NormBadness();
	    if (! lo_div->full() )
	    {
		// save clone of advancing winner
		PMOL winner_clone = new Molecule(*advancing);
		lo_div->push_back(winner_clone);
	    }
	    // advance as far as possible
	    advancing->Evolve(rp.dist_trials, rp.tri_trials, rp.pyr_trials,
		    rp.lookout_prob);
	    int hi_level = advancing->NAtoms();
	    if (rp.trace)
	    {
		advancing->trace.push_back(rv.season);
		advancing->trace.push_back(lo_level);
		advancing->trace.push_back(hi_level);
	    }
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
	    if (eps_gt(descending->NormBadness(), advancing->NormBadness()))
	    {
		*descending = *advancing;
	    }
	    descending->Degenerate(hi_level-lo_level);
	    if (rp.trace)
	    {
		if (lo_level == base_level)
		    descending->trace.clear();
		else
		{
		    descending->trace.push_back(rv.season);
		    descending->trace.push_back(hi_level);
		    descending->trace.push_back(lo_level);
		}
	    }
	    // all set now so we can swap winner and looser
	    (*hi_div)[looser_idx] = advancing;
	    (*lo_div)[winner_idx] = descending;
	    // make sure the original best cluster was much much better than
	    // whatever left in the low division
	    const double spoil_factor = 10.0;
	    if ( advancing_best && eps_gt(lo_div->best()->NormBadness(),
			spoil_factor*adv_bad0) )
	    {
		PMOL lo_looser = lo_div->at(lo_div->find_looser());
		*lo_looser = *advancing;
		for (int nlast = hi_level-1; nlast >= lo_level; --nlast)
		    lo_looser->Pop(nlast);
	    }
	    if (!rp.quiet)
	    {
		cout << rv.season;
		cout << " A " <<
		    lo_level << ' ' << adv_bad0 << ' ' <<
		    hi_level << ' ' << advancing->NormBadness() << "    ";
		cout << " D " <<
		    hi_level << ' ' << desc_bad0 << ' ' <<
		    lo_level << ' ' << descending->NormBadness() << endl;
	    }
	    save_frames_trace(advancing, descending, rp, rv);
	    // no need to finish the season if we found champion
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
	cout << rv.season << " WC " << world_champ->NAtoms() << ' '
	    << world_champ->NormBadness() << endl;
	if (    world_champ->NAtoms() > best_champ.NAtoms() ||
		eps_lt(world_champ->NormBadness(), best_champ.NormBadness())
	   )
	{
	    best_champ = *world_champ;
	    cout << rv.season << " BC " << best_champ.NAtoms() << ' '
		<< best_champ.NormBadness() << endl;
	}
        save_outstru(liga, rp, rv);
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
	exit_code = EXIT_FAILURE;
    }
    else
    {
	cout << "Solution found!!!" << endl << endl;
	exit_code = EXIT_SUCCESS;
    }
    if (rp.trace)
    {
	// needs clean up
	for ( list<int>::iterator ii = best_champ.trace.begin();
		ii != best_champ.trace.end(); )
	{
	    cout << "TR";
	    for (int n = 0; n != 3 && ii != best_champ.trace.end(); ++n)
		cout << ' ' << *(ii++);
	    cout << endl;
	}
	cout << endl;
    }
    BGA::cnt.PrintRunStats();
    // save last frame
    rv.exiting = true;
    save_outstru(liga, rp, rv);
    save_frames(best_champ, rp, rv);
    if (SIGHUP_received)
    {
	exit(exit_code);
    }
    return exit_code;
}
