/***********************************************************************
* Short Title: molecule reconstruction from distance table
*
* Comments: self-tuning competition search for molecular configuration
*
* $Id$
***********************************************************************/

#include <limits>
#include <unistd.h>
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
    // /usage:/;/;/-s/.*/"&\\n"/
    // /cou/;/;/s/^\s*"\(.*\)\\n"/\1/ | '[put! ='/*' | /;/put ='*/'
    cout << 
"usage: " << a.cmd_t << "[-p PAR_FILE] [DISTFILE] [par1=val1 par2=val2...]\n"
"run gizaliga simulation using distances from DISTFILE.  Parameters can\n"
"be set in PAR_FILE or on the command line, which overrides PAR_FILE.\n"
"Options:\n"
"  -p, --parfile=FILE    read parameters from FILE\n"
"  -h, --help            display this message\n"
"  -v, --version         show program version\n"
"IO parameters:\n"
"  distfile=FILE         target distance table\n"
"  outstru=FILE          where to save the best full molecule\n"
"  inistru=FILE          initial structure [empty box]\n"
"  snapshot=FILE         live molecule structure\n"
"  snaprate=int          [100] number of iterations between snapshot updates\n"
"  frames=FILE           save intermediate structures to FILE.liga_round\n"
"  framesrate=int        [100] number of iterations between frame saves\n"
"Liga parameters\n"
"  tol_dd=double         [0.1] distance is not used when dd=|d-d0|>=tol_dd\n"
"  tol_bad=double        [1E-4] target normalized molecule badness\n"
"  seed=int              seed of random number generator\n"
"  ligasize=int          [10] number of teams per division\n"
"  stopgame=double       [1.0] badness threshold to start new round\n"
"  penalty=string        dd penalty function [pow2], fabs, well\n"
"  dist_trials=int       [10] good distance atoms to try\n"
"  tri_trials=int        [20] godd triangle atoms to try\n"
"  pyr_trials=int        [1000] good pyramid atoms to try\n"
;
}

struct RunPar_t
{
    // IO parameters
    string distfile;
    string outstru;
    string inistru;
    string snapshot;
    int snaprate;
    string frames;
    int framesrate;
    // Walk parameters
    double tol_dd;
    double tol_bad;
    int seed;
    int ligasize;
    double stopgame;
    string penalty;
    int dist_trials;
    int tri_trials;
    int pyr_trials;
};

struct RunVar_t
{
    RunVar_t() : liga_round(0), full_liga(false), lastframe(false)
    { }
    int liga_round;
    bool full_liga;
    bool lastframe;
};

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



Molecule process_arguments(RunPar_t& rp, int argc, char *argv[])
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
        print_version();
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
    // assign run parameters
    // distfile
    if (a.args.size())
        a.pars["distfile"] = a.args[0];
    if (!a.ispar("distfile"))
    {
        cerr << "Distance file not defined" << endl;
        exit(EXIT_FAILURE);
    }
    rp.distfile = a.pars["distfile"];
    DistanceTable* dtab;
    try {
        dtab = new DistanceTable(rp.distfile.c_str());
    }
    catch (IOError) {
        exit(EXIT_FAILURE);
    }
    string hashsep(72, '#');
    cout << hashsep << endl;
    cout << "# " << a.cmd_t << ' ' <<
        "$Id$" << endl;
    time_t cur_time = time(NULL);
    cout << "# " << ctime(&cur_time);
    cout << hashsep << endl;
    Molecule mol(*dtab);
    cout << "distfile=" << rp.distfile << endl;
    if (a.ispar("outstru"))
    {
        rp.outstru = a.pars["outstru"];
        cout << "outstru=" << rp.outstru << endl;
    }
    if (a.ispar("inistru"))
    {
        rp.inistru = a.pars["inistru"];
        cout << "inistru=" << rp.inistru << endl;
        try {
            mol.ReadXYZ(rp.inistru.c_str());
        }
        catch (IOError) {
            exit(EXIT_FAILURE);
        }
    }
    if (a.ispar("snapshot"))
    {
        rp.snapshot = a.pars["snapshot"];
        cout << "snapshot=" << rp.snapshot << endl;
        rp.snaprate = a.GetPar<int>("snaprate", 100);
        cout << "snaprate=" << rp.snaprate << endl;
    }
    if (a.ispar("frames"))
    {
        rp.frames = a.pars["frames"];
        cout << "frames=" << rp.frames << endl;
        rp.framesrate = a.GetPar<int>("framesrate", 100);
        cout << "framesrate=" << rp.framesrate << endl;
    }
    // Walk parameters
    rp.tol_dd = a.GetPar<double>("tol_dd", 0.1);
    cout << "tol_dd=" << rp.tol_dd << endl;
    mol.tol_dd = rp.tol_dd;
    rp.tol_bad = a.GetPar<double>("tol_bad", 1.0e-4);
    cout << "tol_bad=" << rp.tol_bad << endl;
    rp.seed = a.GetPar<int>("seed", 0);
    if (rp.seed)
    {
        gsl_rng_set(BGA::rng, rp.seed);
        cout << "seed=" << rp.seed << endl;
    }
    rp.ligasize = a.GetPar<int>("ligasize", 10);
    cout << "ligasize=" << rp.ligasize << endl;
    rp.stopgame = a.GetPar<double>("stopgame", 1.0);
    cout << "stopgame=" << rp.stopgame << endl;
    rp.penalty = a.GetPar<string>("penalty", "pow2");
    if (rp.penalty == "pow2")
        mol.penalty = BGA::pow2;
    else if (rp.penalty == "well")
        mol.penalty = BGA::well;
    else if (rp.penalty == "fabs")
        mol.penalty = fabs;
    else
    {
        cerr << "Invalid value of penalty parameter" << endl;
        exit(EXIT_FAILURE);
    }
    cout << "penalty=" << rp.penalty << endl;
    rp.dist_trials = a.GetPar("dist_trials", 10);
    cout << "dist_trials=" << rp.dist_trials << endl;
    rp.tri_trials = a.GetPar("tri_trials", 20);
    cout << "tri_trials=" << rp.tri_trials << endl;
    rp.pyr_trials = a.GetPar("pyr_trials", 1000);
    cout << "pyr_trials=" << rp.pyr_trials << endl;
    cout << hashsep << endl << endl;
    return mol;
}

void save_snapshot(Molecule& mol, RunPar_t& rp)
{
//  numeric_limits<double> double_info;
    static int cnt = 0;
    static double bestMNB = numeric_limits<double>().max();
    if (rp.snapshot.size() == 0 || rp.snaprate == 0 || ++cnt < rp.snaprate)
        return;
    if ( mol.NormBadness() < bestMNB)
    {
        bestMNB = mol.NormBadness();
        mol.WriteAtomEye(rp.snapshot.c_str());
        cnt = 0;
    }
}

void save_frames(Molecule& mol, RunPar_t& rp, RunVar_t& rv)
{
//  numeric_limits<double> double_info;
    static int cnt = 0;
    if (  rp.frames.size() == 0 || rp.framesrate == 0 ||
            (++cnt < rp.framesrate && !rv.lastframe) )
        return;
    ostringstream oss;
    oss << rp.frames << "." << rv.liga_round;
    mol.WriteAtomEye(oss.str().c_str());
    cnt = 0;
}

int main(int argc, char *argv[])
{
    // process arguments
    RunPar_t rp;
    RunVar_t rv;
    Molecule mol = process_arguments(rp, argc, argv);

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
	<< first_team->Badness() << endl;
    liga[mol.NAtoms()].push_back(first_team);
    // fill lower divisions
    cout << "Filling lower divisions" << endl;
    for (int level = mol.NAtoms()-1; level >= 0; --level)
    {
        PMOL parent_team = liga[level+1].back();
        PMOL lower_team = new Molecule(*parent_team);
        lower_team->Degenerate(1);
	cout << rv.liga_round << "L " << lower_team->NAtoms() << ' '
	    << lower_team->Badness() << endl;
        liga[level].push_back(lower_team);
    }
    cout << "Done" << endl;
    // fill higher divisions
    cout << "Filling lower divisions" << endl;
    mol.evolve_jump = false;
    for (int level = mol.NAtoms()+1; level <= mol.max_NAtoms(); ++level)
    {
        PMOL parent_team = liga[level-1].back();
        PMOL higher_team = new Molecule(*parent_team);
        higher_team->Evolve(rp.dist_trials, rp.tri_trials, rp.pyr_trials);
	cout << rv.liga_round << " I " << higher_team->NAtoms() << ' '
	    << higher_team->Badness() << endl;
        liga[level].push_back(higher_team);
    }
    mol.evolve_jump = true;
    // find the first world champion
    PMOL world_champ = liga.back().best();
    cout << rv.liga_round << " WC " << world_champ->NAtoms() << ' '
	<< world_champ->Badness() << endl;
    // let the game begin
    while ( !(world_champ->Badness() < rp.tol_bad ) )
    {
        ++rv.liga_round;
	typedef vector<Division_t>::iterator VDit;
	int lo_level = 0;
	for (VDit lo_div = liga.begin();
		lo_div < liga.end()-1; ++lo_div, ++lo_level)
	{
	    int winner_idx = lo_div->find_winner();
	    PMOL advancing = lo_div->at(winner_idx);
	    if (! (advancing->Badness() < rp.stopgame) )
		continue;
	    double adv_bad0 = advancing->Badness();
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
	    int looser_idx = hi_div->find_looser();
	    PMOL descending = hi_div->at(looser_idx);
	    double desc_bad0 = descending->Badness();
	    if (! hi_div->full() )
	    {
		// save clone of descending looser
		PMOL looser_clone = new Molecule(*descending);
		hi_div->push_back(looser_clone);
	    }
	    descending->Degenerate(hi_level-lo_level);
	    // all set now so we can swap winner and looser
	    (*hi_div)[looser_idx] = advancing;
	    (*lo_div)[winner_idx] = descending;
	    cout << rv.liga_round;
	    cout << " A " <<
		lo_level << ' ' << adv_bad0 << ' ' <<
		hi_level << ' ' << advancing->Badness() << "    ";
	    cout << " D " <<
		hi_level << ' ' << desc_bad0 << ' ' <<
		lo_level << ' ' << descending->Badness() << endl;
	}
	world_champ = liga.back().best();
	cout << rv.liga_round << " WC " << world_champ->NAtoms() << ' '
	    << world_champ->Badness() << endl;
        save_snapshot(*world_champ, rp);
        save_frames(*world_champ, rp, rv);
    }
    cout << "Solution found!!!" << endl;
    cout << "cnt_penalty_calls = " <<
	BGA::cnt_penalty_calls << endl;
    // save last frame
    rv.lastframe = true;
    save_frames(mol, rp, rv);
    // save final structure
    if (rp.outstru.size() != 0)
        world_champ->WriteAtomEye(rp.outstru.c_str());
    return EXIT_SUCCESS;
}
