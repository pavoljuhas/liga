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
#include "Random.hpp"
#include "TrialDistributor.hpp"
#include "StringUtils.hpp"
#include "Molecule.hpp"
#include "Liga_t.hpp"
#include "AtomFilter_t.hpp"
#include "Exceptions.hpp"
#include "RegisterSVNId.hpp"

using namespace std;
using namespace LIGA;

RegisterSVNId RunPar_t_cpp_id("$Id$");


////////////////////////////////////////////////////////////////////////
// class RunPar_t
////////////////////////////////////////////////////////////////////////

RunPar_t::RunPar_t(int argc, char* argv[])
{
    processArguments(argc, argv);
}

void RunPar_t::processArguments(int argc, char* argv[])
{
    char *short_options =
	"p:hV";
    // parameters and options
    option long_options[] = {
	{"parfile", 1, 0, 'p'},
	{"help", 0, 0, 'h'},
	{"version", 0, 0, 'V'},
	{0, 0, 0, 0}
    };
    args.reset(new ParseArgs(argc, argv, short_options, long_options));
    args->Parse();
    if (args->isopt("h") || argc == 1)
    {
	print_help();
	exit(EXIT_SUCCESS);
    }
    else if (args->isopt("V"))
    {
	cout << version_string() << endl;
	exit(EXIT_SUCCESS);
    }
    if (args->isopt("p"))
    {
	args->ReadPars(args->opts["p"].c_str());
    }
    args->ValidatePars(validpars());
    // assign run parameters
    // distfile
    if (args->args.size())
    {
	args->pars["distfile"] = args->args[0];
    }
    if (args->args.size() > 1)
    {
	ostringstream emsg;
	emsg << argv[0] << ": several DISTFILE arguments.";
	throw ParseArgsError(emsg.str());
    }
    if (!args->ispar("distfile"))
    {
	string emsg = "Distance file not defined.";
	throw ParseArgsError(emsg);
    }
    distfile = args->pars["distfile"];
    DistanceTable dtab;
    ifstream dstfid(distfile.c_str());
    if (!dstfid)
    {
        ostringstream emsg;
        emsg << "Unable to read '" << distfile << "'";
        throw IOError(emsg.str());
    }
    dstfid >> dtab;
    // create empty molecule
    mol.reset(new Molecule(dtab));
    // inistru
    if (args->ispar("inistru"))
    {
	inistru = args->pars["inistru"];
	mol->ReadXYZ(inistru.c_str());
    }
    // outstru
    if (args->ispar("outstru"))
    {
	outstru = args->pars["outstru"];
    }
    // outfmt
    outfmt = args->GetPar<string>("outfmt", "xyz");
    if (outfmt == "xyz")
    {
        Molecule::OutFmtXYZ();
    }
    else if (outfmt == "atomeye")
    {
        Molecule::OutFmtAtomEye();
    }
    else
    {
        const char* emsg = "Invalid value of outfmt parameter.";
        throw ParseArgsError(emsg);
    }
    // saverate, saveall
    if (args->ispar("outstru"))
    {
	// saverate
	saverate = args->GetPar<int>("saverate", 10);
	// saveall
	saveall = args->GetPar<bool>("saveall", false);
    }
    // frames, framestrace, framesrate
    if (args->ispar("frames"))
    {
	// frames
	frames = args->pars["frames"];
	// framestrace
	if (args->ispar("framestrace"))
	{
	    framesrate = 0;
	    vector<long> stp = args->GetParVec<long>("framestrace");
	    if (stp.size() % 3)
	    {
		string emsg = "framestrace must have 3n entries.";
		throw ParseArgsError(emsg);
	    }
	    // check if seasons in stp are ordered
	    for (size_t i = 3; i < stp.size(); i += 3)
	    {
		if (stp[i] < stp[i-3])
		{
		    string emsg = "framestrace seasons must be sorted.";
		    throw ParseArgsError(emsg);
		}
	    }
	    TraceId_t tid;
	    for (size_t i = 0; i < stp.size(); i += 3)
	    {
		tid.season = stp[i];
		tid.level = stp[i+1];
		tid.mol_id = stp[i+2];
		framestrace.push_back(tid);
	    }
	}
	// framesrate - parse only if framestrace is empty
	if (framestrace.size() == 0)
	{
	    framesrate = args->GetPar<int>("framesrate", 0);
	}
    }
    // trace
    trace = args->GetPar<bool>("trace", false);
    // verbose
    if (args->ispar("verbose"))
    {
	verbose = args->GetParVec<string>("verbose");
	vector<string>::iterator w;
	for (w = verbose.begin(); w != verbose.end(); ++w)
	{
	    if (!Liga_t::isVerboseFlag(*w))
	    {
		ostringstream emsg;
		emsg << "verbose flag must be one of (" <<
		    joined_verbose_flags() << ")";
		throw ParseArgsError(emsg.str());
	    }
	}
    }
    // liga parameters
    // ndim
    ndim = args->GetPar<size_t>("ndim", 3);
    if (ndim < 1 || ndim > 3)
    {
	string emsg = "ndim value must be 1, 2 or 3.";
	throw ParseArgsError(emsg);
    }
    // crystal
    crystal = args->GetPar<bool>("crystal", false);
    // latpar
    latpar.resize(6);
    latpar[0] = latpar[1] = latpar[2] = 1.0;
    latpar[3] = latpar[4] = latpar[5] = 90.0;
    if (args->ispar("latpar"))
    {
        latpar = args->GetParVec<double>("latpar");
        if (latpar.size() != 6)
        {
            char* emsg = "latpar must define 6 lattice parameters";
            throw ParseArgsError(emsg);
        }
    }
    // distreuse
    distreuse = args->GetPar<bool>("distreuse", false);
    Molecule::distreuse = distreuse;
    // tol_bad
    tol_bad = args->GetPar<double>("tol_bad", 1.0e-4);
    Molecule::tol_nbad = tol_bad;
    // natoms must be set after distreuse
    if (args->ispar("natoms"))
    {
	natoms = args->GetPar<int>("natoms");
	mol->setMaxNAtoms(natoms);
    }
    natoms = mol->maxNAtoms();
    // fixed_atoms must be set after inistru
    if (args->ispar("fixed_atoms"))
    {
	fixed_atoms = args->ExpandRangePar("fixed_atoms");
	sort(fixed_atoms.begin(), fixed_atoms.end());
	vector<int>::iterator last;
	last = unique(fixed_atoms.begin(), fixed_atoms.end());
	fixed_atoms.erase(last, fixed_atoms.end());
	vector<int>::const_iterator ii;
	try {
	    for (ii = fixed_atoms.begin(); ii != fixed_atoms.end(); ++ii)
	    {
		mol->Fix(*ii - 1);
	    }
	}
	catch (range_error) {
	    ostringstream emsg;
	    emsg << "fixed_atoms - invalid index " << *ii << '.';
	    throw ParseArgsError(emsg.str());
	}
    }
    base_level = mol->NFixed();
    // seed_clusters
    if (args->ispar("seed_clusters"))
    {
	vector<int> scs = args->GetParVec<int>("seed_clusters");
	if (scs.size() % 3)
	{
	    string emsg = "seed_clusters must have 3n entries.";
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
	    if (scid.level <= base_level || scid.level > mol->maxNAtoms())
	    {
		ostringstream emsg;
		emsg << "seed_clusters - invalid level " << scid.level << '.';
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
    // maxcputime
    maxcputime = args->GetPar<double>("maxcputime", 0.0);
    // rngseed
    rngseed = args->GetPar<int>("rngseed", 0);
    if (rngseed)
    {
	randomSeed(rngseed);
    }
    // evolve_frac
    evolve_frac = args->GetPar<double>("evolve_frac", 0.1);
    Molecule::evolve_frac = evolve_frac;
    // evolve_relax
    evolve_relax = args->GetPar<bool>("evolve_relax", false);
    Molecule::evolve_relax = evolve_relax;
    // degenerate_relax
    degenerate_relax = args->GetPar<bool>("degenerate_relax", false);
    Molecule::degenerate_relax = degenerate_relax;
    // ligasize
    ligasize = args->GetPar<int>("ligasize", 10);
    // stopgame
    stopgame = args->GetPar<double>("stopgame", 0.0025);
    // seasontrials
    seasontrials = int( args->GetPar<double>("seasontrials", 16384.0) );
    // trials_sharing
    trials_sharing = args->GetPar<string>("trials_sharing", "equal");
    if (!TrialDistributor::isType(trials_sharing))
    {
	ostringstream emsg;
	emsg << "trials_sharing must be one of (";
	emsg << join(", ", TrialDistributor::getTypes()) << ").";
	throw ParseArgsError(emsg.str());
    }
    // lookout_prob
    lookout_prob = args->GetPar("lookout_prob", 0.0);
    Molecule::lookout_prob = lookout_prob;
    // bangle_range
    if (args->ispar("bangle_range"))
    {
	bangle_range = args->GetParVec<double>("bangle_range");
	if (bangle_range.size() != 2 && bangle_range.size() != 3)
	{
	    string emsg = "bangle_range must have 2 or 3 arguments.";
	    throw ParseArgsError(emsg);
	}
	double max_blen = bangle_range[0];
        if (bangle_range.size() == 2)   bangle_range.push_back(DOUBLE_MAX);
	BondAngleFilter_t* pbaf = new BondAngleFilter_t(max_blen);
        pbaf->setBondAngleRange(bangle_range[1], bangle_range[2]);
        Molecule::atom_filters.push_back(pbaf);
    }
    // max_dist
    if (args->ispar("max_dist"))
    {
	max_dist = args->GetPar<double>("max_dist");
	LoneAtomFilter_t* plaf = new LoneAtomFilter_t(max_dist);
        Molecule::atom_filters.push_back(plaf);
    }
    // done
    print_pars();
}

void RunPar_t::print_help()
{
    // /usage:/;/;/-s/.*/"&\\n"/
    // /cou/;/;/s/^\s*"\(.*\)\\n"/\1/ | '[put! ='/*' | /;/put ='*/'
    string cmd_t = args->cmd_t;
    cout << 
"usage: " << cmd_t << " [-p PAR_FILE] [DISTFILE] [par1=val1 par2=val2...]\n"
"run gizaliga simulation using distances from DISTFILE.  Parameters can\n"
"be set in PAR_FILE or on the command line, which overrides PAR_FILE.\n"
"Options:\n"
"  -p, --parfile=FILE    read parameters from FILE\n"
"  -h, --help            display this message\n"
"  -V, --version         show program version\n"
"IO parameters:\n"
"  distfile=FILE         target distance table\n"
"  inistru=FILE          [empty] initial structure in Cartesian coordinates\n"
"  outstru=FILE          where to save the best full molecule\n"
"  outfmt=string         [xyz], atomeye - outstru file format\n"
"  saverate=int          [10] minimum iterations between outstru updates\n"
"  saveall=bool          [false] save best molecules from all levels\n"
"  frames=FILE           save intermediate structures to FILE.season\n"
"  framesrate=int        [0] number of iterations between frame saves\n"
"  framestrace=array     [] triplets of (season, level, id)\n"
"  trace=bool            [false] keep and show trace of the best structure\n"
"  verbose=array         [ad,wc,bc] output flags from (" <<
	joined_verbose_flags() << ")\n" <<
"Liga parameters:\n"
"  ndim={1,2,3}          [3] search in n-dimensional space\n"
"  crystal=bool          [false] assume periodic crystal structure\n"
"  latpar=array          [1,1,1,90,90,90] crystal lattice parameters\n"
"  distreuse=bool        [false] keep used distances in distance table\n"
"  tol_bad=double        [1E-4] target normalized molecule badness\n"
"  natoms=int            use with loose distfiles or for set distreuse\n"
"  fixed_atoms=ranges    [] indices of fixed atoms in inistru (start at 1)\n"
"  seed_clusters=array   [] triplets of (level, number, trials)\n"
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
    oss << quote << "gizaliga " << RegisterSVNId::lastId() << '\n' <<
	quote << "compiler version " << __VERSION__;
    return oss.str();
}

void RunPar_t::print_pars()
{
    // print out all run parameters
    // intro messages
    string hashsep(72, '#');
    cout << hashsep << '\n';
    cout << "# " << args->cmd_t << '\n';
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
    if (args->ispar("inistru"))
    {
	cout << "inistru=" << inistru << '\n';
    }
    // outstru
    if (args->ispar("outstru"))
    {
        cout << "outstru=" << outstru << '\n';
    }
    // outfmt
    if (args->ispar("outstru") || args->ispar("frames"))
    {
	cout << "outfmt=" << outfmt << '\n';
    }
    // saverate, saveall
    if (args->ispar("outstru"))
    {
        cout << "saverate=" << saverate << '\n';
	cout << "saveall=" << saveall << '\n';
    }
    // frames, framestrace, framesrate
    if (args->ispar("frames"))
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
		cout << ' ' << ii->mol_id;
	    }
	    cout << '\n';
	}
	// framesrate - print only if framestrace is empty
	else
	{
	    cout << "framesrate=" << framesrate << '\n';
	}
    }
    // trace
    cout << "trace=" << trace << '\n';
    // verbose
    {
	cout << "verbose=" << join(",", verbose) << '\n';
    }
    // liga parameters
    // ndim, distreuse, tol_bad 
    cout << "ndim=" << ndim << '\n';
    cout << "distreuse=" << distreuse << '\n';
    cout << "tol_bad=" << tol_bad << '\n';
    // natoms
    cout << "natoms=" << natoms << '\n';
    // fixed_atoms
    if (args->ispar("fixed_atoms") && !fixed_atoms.empty())
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
    if (args->ispar("bangle_range"))
    {
	cout << "bangle_range=" << bangle_range[0];
	for (size_t i = 1; i < bangle_range.size(); ++i)
	{
	    cout << ',' << bangle_range[i];
	}
	cout << '\n';
    }
    // max_dist
    if (args->ispar("max_dist"))
    {
	cout << "max_dist=" << max_dist << '\n';
    }
    // finish done
    cout << hashsep << '\n' << endl;
}

const list<string>& RunPar_t::validpars() const
{
    char *parnames[] = {
        "distfile",
        "inistru",
        "outstru",
        "outfmt",
        "saverate",
        "saveall",
        "frames",
        "framesrate",
        "framestrace",
        "verbose",
        "ndim",
        "distreuse",
        "tol_bad",
        "natoms",
        "fixed_atoms",
        "seed_clusters",
        "maxcputime",
        "rngseed",
        "evolve_frac",
        "evolve_relax",
        "degenerate_relax",
        "ligasize",
        "stopgame",
        "seasontrials",
	"trace",
        "trials_sharing",
        "lookout_prob",
        "bangle_range",
        "max_dist",
    };
    size_t parlen = sizeof(parnames)/sizeof(char*);
    static list<string> parlist(parnames, parnames + parlen);
    return parlist;
}

const string& RunPar_t::joined_verbose_flags() const
{
    static string jvf;
    if (jvf.empty())
    {
        jvf = join(", ", Liga_t::verbose_flags);
    }
    return jvf;
}

// End of file
