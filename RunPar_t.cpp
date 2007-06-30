/***********************************************************************
* Short Title: run parameters for gizaliga application
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include "RunPar_t.hpp"
#include "TrialDistributor.hpp"
#include "StringUtils.hpp"

using namespace std;

RegisterSVNId RunPar_t_cpp_id("$Id$");


////////////////////////////////////////////////////////////////////////
// class RunPar_t
////////////////////////////////////////////////////////////////////////

RunPar_t::RunPar_t()
{
    fill_validpars();
    fill_verbose();
    crystal = false;
    latpar.resize(6);
    latpar[0] = latpar[1] = latpar[2] = 1.0;
    latpar[3] = latpar[4] = latpar[5] = 90.0;
}

void RunPar_t::processArguments(int argc, char *argv[])
{
    char *short_options =
	"p:qthV";
    // parameters and options
    option long_options[] = {
	{"parfile", 1, 0, 'p'},
	{"trace", 0, 0, 't'},
	{"help", 0, 0, 'h'},
	{"version", 0, 0, 'V'},
	{0, 0, 0, 0}
    };
    ParseArgs a(argc, argv, short_options, long_options);
    a.Parse();
    if (a.isopt("h") || argc == 1)
    {
	print_help(a);
	exit(EXIT_SUCCESS);
    }
    else if (a.isopt("V"))
    {
	cout << version_string() << endl;
	exit(EXIT_SUCCESS);
    }
    if (a.isopt("p"))
    {
	a.ReadPars(a.opts["p"].c_str());
    }
    trace = a.isopt("t");
    a.ValidatePars(validpars);
    // assign run parameters
    // distfile
    if (a.args.size())
	a.pars["distfile"] = a.args[0];
    if (a.args.size() > 1)
    {
	ostringstream emsg;
	emsg << argv[0] << ": several DISTFILE arguments";
	throw ParseArgsError(emsg.str());
    }
    if (!a.ispar("distfile"))
    {
	string msg = "Distance file not defined";
	throw ParseArgsError(msg);
    }
    distfile = a.pars["distfile"];
    DistanceTable* dtab;
    dtab = new DistanceTable(distfile.c_str());
    // create empty molecule
    mol = Molecule(*dtab);
    // inistru
    if (a.ispar("inistru"))
    {
	inistru = a.pars["inistru"];
	mol.ReadXYZ(inistru.c_str());
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
	else throw ParseArgsError("Invalid value of outfmt parameter");
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
		string emsg = "framestrace must have 3n entries";
		throw ParseArgsError(emsg);
	    }
	    // check if seasons in stp are ordered
	    for (size_t i = 3; i < stp.size(); i += 3)
	    {
		if (stp[i] < stp[i-3])
		{
		    string emsg = "framestrace seasons must be sorted";
		    throw ParseArgsError(emsg);
		}
	    }
	    TraceId_t tid;
	    for (size_t i = 0; i < stp.size(); i += 3)
	    {
		tid.season = stp[i];
		tid.level = stp[i+1];
		tid.index = stp[i+2];
		framestrace.push_back(tid);
	    }
	}
	// framesrate - parse only if framestrace is empty
	if (framestrace.size() == 0)
	{
	    framesrate = a.GetPar<int>("framesrate", 0);
	}
    }
    // verbose
    if (a.ispar("verbose"))
    {
	verbose = false;
	vector<string> flagwords = a.GetParVec<string>("verbose");
	vector<string>::iterator w, vbsflag;
	for (w = flagwords.begin(); w != flagwords.end(); ++w)
	{
	    vbsflag = find(verbose_flag.begin(), verbose_flag.end(), *w);
	    if (vbsflag == verbose_flag.end())
	    {
		ostringstream emsg;
		emsg << "verbose flag must be one of (" <<
		    join(", ", verbose_flag) << ")";
		throw ParseArgsError(emsg.str());
	    }
	    int vbsindex = vbsflag - verbose_flag.begin();
	    verbose[vbsindex] = true;
	    if (vbsindex == VerboseFlag::ALL)	verbose = true;
	}
	verbose_mute = !verbose.sum();
    }
    // liga parameters
    // ndim
    ndim = a.GetPar<size_t>("ndim", 3);
    if (ndim < 1 || ndim > 3)
    {
	string emsg = "ndim value must be 1, 2 or 3";
	throw ParseArgsError(emsg);
    }
    // crystal
    crystal = a.GetPar<bool>("crystal", false);
    // latpar
    if (a.ispar("latpar"))
    {
        latpar = a.GetParVec<double>("latpar");
        if (latpar.size() != 6)
        {
            char* emsg = "latpar must define 6 lattice parameters";
            throw ParseArgsError(emsg);
        }
    }
    // tol_dd
    tol_dd = a.GetPar<double>("tol_dd", 0.1);
    mol.tol_dd = tol_dd;
    // tol_bad
    tol_bad = a.GetPar<double>("tol_bad", 1.0e-4);
    mol.tol_nbad = tol_bad;
    // natoms must be set after tol_dd
    if (a.ispar("natoms"))
    {
	natoms = a.GetPar<int>("natoms");
	mol.setMaxNAtoms(natoms);
    }
    natoms = mol.maxNAtoms();
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
	    ostringstream emsg;
	    emsg << "fixed_atoms - invalid index " << *ii;
	    throw ParseArgsError(emsg.str());
	}
    }
    base_level = mol.NFixed();
    // seed_clusters
    if (a.ispar("seed_clusters"))
    {
	vector<int> scs = a.GetParVec<int>("seed_clusters");
	if (scs.size() % 3)
	{
	    string emsg = "seed_clusters must have 3n entries";
	    throw ParseArgsError(emsg);
	}
	map<int,SeedClusterInfo> lvsc;
	SeedClusterInfo scid;
	// make sure seed_clusters are sorted and level-unique
	for (size_t i = 0; i < scs.size(); i += 3)
	{
	    scid.level = scs[i];
	    scid.number = scs[i+1];
	    scid.trials = scs[i+2];
	    if (scid.level <= base_level || scid.level > mol.maxNAtoms())
	    {
		ostringstream emsg;
		emsg << "seed_clusters - invalid level " << scid.level;
		throw ParseArgsError(emsg.str());
	    }
	    lvsc[scid.level] = scid;
	}
	for (   map<int,SeedClusterInfo>::iterator ii = lvsc.begin();
		ii != lvsc.end(); ++ii )
	{
	    seed_clusters.push_back(ii->second);
	}
    }
    // centersize
    centersize = a.GetPar<int>("centersize", 0);
    mol.center_size = centersize;
    // maxcputime
    maxcputime = a.GetPar<double>("maxcputime", 0.0);
    // rngseed
    rngseed = a.GetPar<int>("rngseed", 0);
    if (rngseed)
    {
	gsl_rng_set(BGA::rng, rngseed);
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
    // seasontrials
    seasontrials = int( a.GetPar<double>("seasontrials", 16384.0) );
    // trials_sharing
    trials_sharing = a.GetPar<string>("trials_sharing", "equal");
    if (!TrialDistributor::isType(trials_sharing))
    {
	ostringstream emsg;
	emsg << "trials_sharing must be one of (";
	emsg << join(", ", TrialDistributor::getTypes()) << ").";
	throw ParseArgsError(emsg.str());
    }
    // lookout_prob
    lookout_prob = a.GetPar("lookout_prob", 0.0);
    mol.lookout_prob = lookout_prob;
    // bangle_range
    if (a.ispar("bangle_range"))
    {
	bangle_range = a.GetParVec<double>("bangle_range");
	if (bangle_range.size() != 2 && bangle_range.size() != 3)
	{
	    string emsg = "bangle_range must have 2 or 3 arguments";
	    throw ParseArgsError(emsg);
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
    // done
    print_pars(a);
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
"  -t, --trace           keep and show frame trace of the champion\n"
"  -h, --help            display this message\n"
"  -V, --version         show program version\n"
"IO parameters:\n"
"  distfile=FILE         target distance table\n"
"  inistru=FILE          [empty] initial structure in Cartesian coordinates\n"
"  outstru=FILE          where to save the best full molecule\n"
"  outfmt=string         [xyz], atomeye - outstru file format\n"
"  saverate=int          [10] minimum iterations between outstru updates\n"
"  saveall=bool          [false] save best molecules from all divisions\n"
"  frames=FILE           save intermediate structures to FILE.season\n"
"  framesrate=int        [0] number of iterations between frame saves\n"
"  framestrace=array     [] triplets of (season, initial, final level)\n"
"  verbose=array         [AD,WC,BC] output flags from (" <<
	join(",", verbose_flag) << ")\n" <<
"Liga parameters:\n"
"  ndim={1,2,3}          [3] search in n-dimensional space\n"
"  crystal=bool          [false] assume periodic crystal structure\n"
"  latpar=array          [1,1,1,90,90,90] crystal lattice parameters\n"
"  tol_dd=double         [0.1] distance is not used when dd=|d-d0|>=tol_dd\n"
"  tol_bad=double        [1E-4] target normalized molecule badness\n"
"  natoms=int            use with loose distfiles or when tol_dd==0\n"
"  fixed_atoms=ranges    [] indices of fixed atoms in inistru (start at 1)\n"
"  seed_clusters=array   [] triplets of (level, number, trials)\n"
"  centersize=int        [0] shift smaller molecules to the origin\n"
"  maxcputime=double     [0] when set, maximum CPU time in seconds\n"
"  rngseed=int           seed of random number generator\n"
"  evolve_frac=double    [0.1] fraction of tol_bad threshold of tested atoms\n"
"  evolve_relax=bool     [false] relax the worst atom after addition\n"
"  degenerate_relax=bool [false] relax the worst atom after removal\n"
"  ligasize=int          [10] number of teams per division\n"
"  stopgame=double       [0.0025] skip division when winner is worse\n"
"  seasontrials=int      [16384] number of atom placements in one season\n"
"  trials_sharing=string [equal] sharing method from (" <<
	join(",", TrialDistributor::getTypes()) << ")\n" <<
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
    oss << quote << "gizaliga " << RegisterSVNId::last_id << '\n' <<
	quote << "compiler version " << __VERSION__;
    return oss.str();
}

void RunPar_t::print_pars(ParseArgs& a)
{
    // print out all run parameters
    // intro messages
    string hashsep(72, '#');
    cout << hashsep << '\n';
    cout << "# " << a.cmd_t << '\n';
    cout << version_string("# ") << '\n';
    char hostname[255];
    gethostname(hostname, 255);
    cout << "# " << hostname << '\n';
    time_t cur_time = time(NULL);
    cout << "# " << ctime(&cur_time);
    cout << hashsep << '\n';
    // distfile
    cout << "distfile=" << distfile << '\n';
    // inistru
    if (a.ispar("inistru"))
    {
	cout << "inistru=" << inistru << '\n';
    }
    // outstru, outfmt, saverate, saveall
    if (a.ispar("outstru"))
    {
        cout << "outstru=" << outstru << '\n';
	cout << "outfmt=" << outfmt << '\n';
        cout << "saverate=" << saverate << '\n';
	cout << "saveall=" << saveall << '\n';
    }
    // frames, framestrace, framesrate
    if (a.ispar("frames"))
    {
	cout << "frames=" << frames << '\n';
	// framestrace
	if (framestrace.size() != 0)
	{
	    cout << "framestrace=";
	    for (deque<TraceId_t>::iterator ii = framestrace.begin();
		    ii != framestrace.end(); ++ii)
	    {
		cout << " \\" << '\n';
		cout << "    " << ii->season;
		cout << ' ' << ii->level;
		cout << ' ' << ii->index;
	    }
	    cout << '\n';
	}
	// framesrate - print only if framestrace is empty
	else
	{
	    cout << "framesrate=" << framesrate << '\n';
	}
    }
    // verbose
    if (a.ispar("verbose"))
    {
	list<string> flagwords;
	for (size_t i = 0; i < verbose.size(); ++i)
	{
	    if (verbose[i]) flagwords.push_back(verbose_flag[i]);
	}
	cout << "verbose=" << join(",", flagwords) << '\n';
    }
    // liga parameters
    // ndim, tol_dd, tol_bad 
    cout << "ndim=" << ndim << '\n';
    cout << "tol_dd=" << tol_dd << '\n';
    cout << "tol_bad=" << tol_bad << '\n';
    // natoms
    if (a.ispar("natoms"))
    {
	cout << "natoms=" << natoms << '\n';
    }
    // fixed_atoms
    if (a.ispar("fixed_atoms") && !fixed_atoms.empty())
    {
	cout << "fixed_atoms=";
	for (int start = 0; start < int(fixed_atoms.size()); )
	{
	    int range = 0;
	    for ( ; start + range < int(fixed_atoms.size()) &&
		    fixed_atoms[start]+range == fixed_atoms[start+range];
		    ++range )
	    { }
	    if (start > 0)  cout << ',';
	    cout << fixed_atoms[start];
	    if (range > 1)  cout << ".." << fixed_atoms[start+range-1];
	    start += range;
	}
	cout << '\n';
    }
    // seed_clusters
    if (seed_clusters.size() != 0)
    {
	cout << "seed_clusters=";
	for (   vector<SeedClusterInfo>::iterator ii = seed_clusters.begin();
		ii != seed_clusters.end(); ++ii)
	{
	    cout << " \\" << '\n';
	    cout << "    " << ii->level;
	    cout << ' ' << ii->number;
	    cout << ' ' << ii->trials;
	}
	cout << '\n';
    }
    // centersize
    if (a.ispar("centersize"))
    {
	cout << "centersize=" << centersize << '\n';
    }
    // maxcputime
    if (maxcputime > 0.0)
    {
	cout << "maxcputime=" << maxcputime << '\n';
    }
    // rngseed
    if (rngseed)
    {
	cout << "rngseed=" << rngseed << '\n';
    }
    // evolve_frac, evolve_relax, degenerate_relax
    cout << "evolve_frac=" << evolve_frac << '\n';
    cout << "evolve_relax=" << evolve_relax << '\n';
    cout << "degenerate_relax=" << degenerate_relax << '\n';
    // ligasize, stopgame, seasontrials, trials_sharing, lookout_prob
    cout << "ligasize=" << ligasize << '\n';
    cout << "stopgame=" << stopgame << '\n';
    cout << "seasontrials=" << seasontrials << '\n';
    cout << "trials_sharing=" << trials_sharing << '\n';
    cout << "lookout_prob=" << lookout_prob << '\n';
    // constraints
    // bangle_range
    if (a.ispar("bangle_range"))
    {
	cout << "bangle_range=" << bangle_range[0];
	for (size_t i = 1; i < bangle_range.size(); ++i)
	{
	    cout << ',' << bangle_range[i];
	}
	cout << '\n';
    }
    // max_dist
    if (a.ispar("max_dist"))
    {
	cout << "max_dist=" << max_dist << '\n';
    }
    // finish done
    cout << hashsep << '\n' << endl;
}

void RunPar_t::fill_validpars()
{
    char *pnames[] = {
	"distfile", "inistru",
	"outstru", "outfmt", "saverate", "saveall",
	"frames", "framesrate", "framestrace", "verbose", "ndim",
	"tol_dd", "tol_bad", "natoms", "fixed_atoms", "seed_clusters",
	"centersize", "maxcputime", "rngseed",
	"evolve_frac", "evolve_relax", "degenerate_relax",
	"ligasize", "stopgame",
	"seasontrials", "trials_sharing", "lookout_prob",
	"bangle_range", "max_dist" };
    validpars.insert(validpars.end(),
	    pnames, pnames+sizeof(pnames)/sizeof(char*));
}

void RunPar_t::fill_verbose()
{
    using namespace VerboseFlag;
    // default verbose value
    verbose.resize(VERBOSE_SIZE, false);
    verbose[AD] = true;
    verbose[WC] = true;
    verbose[BC] = true;
    verbose_mute = false;
    // fill verbose_flag array
    verbose_flag.resize(VERBOSE_SIZE);
    verbose_flag[AD] = "ad";
    verbose_flag[WC] = "wc";
    verbose_flag[BC] = "bc";
    verbose_flag[AV] = "av";
    verbose_flag[TS] = "ts";
    verbose_flag[ALL] = "all";
}


// End of file
