/***********************************************************************
* Short Title: run parameters for gizaliga application
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <cstdlib>
#include <fstream>
#include "RunPar_t.hpp"
#include "Random.hpp"
#include "TrialDistributor.hpp"
#include "StringUtils.hpp"
#include "Molecule.hpp"
#include "Crystal.hpp"
#include "Lattice.hpp"
#include "Liga_t.hpp"
#include "AtomFilter_t.hpp"
#include "LigaUtils.hpp"
#include "Exceptions.hpp"
#include "Version.hpp"

using namespace std;
using namespace NS_LIGA;

////////////////////////////////////////////////////////////////////////
// class RunPar_t
////////////////////////////////////////////////////////////////////////

RunPar_t::RunPar_t(int argc, char* argv[])
{
    processArguments(argc, argv);
    this->mol->CheckIntegrity();
}

void RunPar_t::processArguments(int argc, char* argv[])
{
    const char* short_options = "p:hV";
    // parameters and options
    option long_options[] = {
	{"parfile", 1, 0, 'p'},
	{"help", 0, 0, 'h'},
	{"version", 0, 0, 'V'},
	{0, 0, 0, 0}
    };
    args.reset(new ParseArgs(argc, argv, short_options, long_options));
    // define aliases for renamed parameters
    args->defParameterAlias("tol_bad", "tolcost");
    args->defParameterAlias("evolve_frac", "promotefrac");
    args->defParameterAlias("evolve_relax", "promoterelax");
    args->defParameterAlias("degenerate_relax", "demoterelax");
    args->defParameterAlias("trials_sharing", "trialsharing");
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
    // figure out if we have Molecule or Crystal
    // crystal
    crystal = args->GetPar<bool>("crystal", false);
    // create empty molecule
    if (this->crystal)
    {
        Crystal crst;
        this->mol.reset(new Crystal(crst));
    }
    else
    {
        mol.reset(new Molecule);
    }
    // assign distance table
    mol->setDistanceTable(dtab);
    // inistru
    if (args->ispar("inistru"))
    {
	this->inistru = args->pars["inistru"];
	this->mol->ReadFile(inistru);
    }
    // latpar
    // must be processed after reading inistru, which may define lattice
    latpar.resize(6);
    if (this->crystal)
    {
        Crystal& crst = dynamic_cast<Crystal&>(*mol);
        const Lattice& L = crst.getLattice();
        latpar[0] = L.a();
        latpar[1] = L.b();
        latpar[2] = L.c();
        latpar[3] = L.alpha();
        latpar[4] = L.beta();
        latpar[5] = L.gamma();
    }
    if (args->ispar("latpar"))
    {
        if (!this->crystal)
        {
            const char* emsg = "latpar has no sense when crystal=false.";
            throw ParseArgsError(emsg);
        }
        latpar = args->GetParVec<double>("latpar");
        if (latpar.size() != 6)
        {
            const char* emsg = "latpar must define 6 lattice parameters.";
            throw ParseArgsError(emsg);
        }
        Lattice L(latpar[0], latpar[1], latpar[2],
                latpar[3], latpar[4], latpar[5]);
        Crystal& crst = dynamic_cast<Crystal&>(*mol);
        crst.setLattice(L);
    }
    // rmax
    if (args->ispar("rmax"))
    {
        if (!this->crystal)
        {
            const char* emsg = "rmax has no sense when crystal=false.";
            throw ParseArgsError(emsg);
        }
        this->rmax = args->GetPar<double>("rmax");
        Crystal& crst = dynamic_cast<Crystal&>(*mol);
        crst.setRmax(rmax);
    }
    if (crystal)
    {
        Crystal& crst = dynamic_cast<Crystal&>(*mol);
        this->rmax = crst.getRmax();
    }
    // outstru
    if (args->ispar("outstru"))
    {
	outstru = args->pars["outstru"];
    }
    // outfmt
    outfmt = args->GetPar<string>("outfmt", "rawxyz");
    Molecule::setOutputFormat(outfmt);
    // saverate, saveall
    if (args->ispar("outstru"))
    {
	// saverate
	saverate = args->GetPar<int>("saverate", 0);
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
                tid.mol_natoms = 0;
                tid.mol_norm_badness = 0.0;
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
    verbose = Liga_t::getDefaultVerbose();
    if (args->ispar("verbose"))
    {
        fill(verbose.begin(), verbose.end(), false);
	vector<string> flagwords = args->GetParVec<string>("verbose");
	vector<string>::iterator w;
        try {
            for (w = flagwords.begin(); w != flagwords.end(); ++w)
            {
                Liga_t::setVerboseVector(verbose, *w, true);
            }
        }
        catch (invalid_argument(e)) {
            throw ParseArgsError(e.what());
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
    // distreuse
    distreuse = args->GetPar<bool>("distreuse", this->crystal);
    mol->setDistReuse(distreuse);
    // tolcost
    tolcost = args->GetPar<double>("tolcost", 1.0e-4);
    Molecule::tol_nbad = tolcost;
    // natoms - apply only when formula is not set
    if (args->ispar("natoms") && !args->ispar("formula"))
    {
        this->natoms = args->GetPar<int>("natoms");
        this->formula.clear();
        this->formula.push_back(
                ChemicalFormula::value_type("C", this->natoms));
        this->mol->setChemicalFormula(this->formula);
    }
    // formula must be set after distreuse
    if (args->ispar("formula"))
    {
	this->formula.fromString(args->GetPar<string>("formula"));
	this->mol->setChemicalFormula(this->formula);
    }
    this->formula = mol->getChemicalFormula();
    // radii
    if (args->ispar("radii"))
    {
        this->radii.fromString(args->GetPar<string>("radii"));
        this->mol->setAtomRadiiTable(this->radii);
    }
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
    // maxcputime
    maxcputime = args->GetPar<double>("maxcputime", 0.0);
    // rngseed
    rngseed = args->GetPar<int>("rngseed", 0);
    if (rngseed)
    {
	randomSeed(rngseed);
    }
    // promotefrac
    promotefrac = args->GetPar<double>("promotefrac", 0.1);
    Molecule::promotefrac = promotefrac;
    // promoterelax
    promoterelax = args->GetPar<bool>("promoterelax", false);
    Molecule::promoterelax = promoterelax;
    // demoterelax
    demoterelax = args->GetPar<bool>("demoterelax", false);
    Molecule::demoterelax = demoterelax;
    // ligasize
    ligasize = args->GetPar<int>("ligasize", 10);
    // stopgame
    stopgame = args->GetPar<double>("stopgame", 0.0025);
    // seasontrials
    seasontrials = int( args->GetPar<double>("seasontrials", 16384.0) );
    // trialsharing
    trialsharing = args->GetPar<string>("trialsharing", "success");
    if (!TrialDistributor::isType(trialsharing))
    {
	ostringstream emsg;
	emsg << "trialsharing must be one of (";
	emsg << join(", ", TrialDistributor::getTypes()) << ").";
	throw ParseArgsError(emsg.str());
    }
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
    const string& cmd_t = args->cmd_t;
    cout << 
"usage: " << cmd_t << " [-p PARFILE] [DISTFILE] [par1=val1 par2=val2...]\n"
"run gizaliga simulation using distances from DISTFILE.  Parameters can\n"
"be set in PARFILE or on the command line, which overrides PARFILE.\n"
"Options:\n"
"  -p, --parfile=FILE    read parameters from FILE\n"
"  -h, --help            display this message\n"
"  -V, --version         show program version\n"
"IO parameters:\n"
"  distfile=FILE         target distance table\n"
"  inistru=FILE          [empty] initial structure in Cartesian coordinates\n"
"  outstru=FILE          where to save the best full molecule\n"
"  outfmt=string         [rawxyz], cif, discus,... - outstru file format\n"
"  saverate=int          [0] rate of intermediate saves of outstru\n"
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
"  rmax=double           [dmax] distance cutoff when crystal=true\n"
"  distreuse=bool        [false] keep used distances in distance table\n"
"  tolcost=double        [1E-4] target normalized molecule cost\n"
"  natoms=int            obsolete, equivalent to formula=Cn\n"
"  formula=string        chemical formula, use inistru when not specified\n"
"  radii=string          define atomic radii in (A1:r1, A2:r2,...) format\n"
"  fixed_atoms=ranges    [] indices of fixed atoms in inistru (start at 1)\n"
"  maxcputime=double     [0] when set, maximum CPU time in seconds\n"
"  rngseed=int           seed of random number generator\n"
"  promotefrac=double    [0.1] fraction of tolcost threshold of tested atoms\n"
"  promoterelax=bool     [false] relax the worst atom after addition\n"
"  demoterelax=bool      [false] relax the worst atom after removal\n"
"  ligasize=int          [10] number of teams per division\n"
"  stopgame=double       [0.0025] skip division when winner is worse\n"
"  seasontrials=int      [16384] number of atom placements in one season\n"
"  trialsharing=string   [success] sharing method from (" <<
	join(",", TrialDistributor::getTypes()) << ")\n" <<
"Constrains (applied only when set):\n"
"  bangle_range=array    (max_blen, low[, high]) bond angle constraint\n"
"  max_dist=double       distance limit for rejecting lone atoms\n"
;
}

string RunPar_t::version_string(string quote)
{
    using namespace std;
    ostringstream oss;
    oss << quote << "gizaliga " << NS_VERSION::getId() << '\n' <<
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
	cout << "inistru=" << this->inistru << '\n';
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
	vector<string> flagwords;
        // scan flags only up to ALL
        using namespace NS_LIGA_VERBOSE_FLAG;
	for (size_t i = 0; i < ALL; ++i)
	{
	    if (verbose[i]) flagwords.push_back(Liga_t::verbose_flags[i]);
	}
	cout << "verbose=" << join(",", flagwords) << '\n';
    }
    // liga parameters
    // ndim
    cout << "ndim=" << ndim << '\n';
    // crystal, latpar, rmax
    cout << "crystal=" << this->crystal << '\n';
    if (this->crystal)
    {
        cout << "latpar=" <<
            this->latpar[0] << ',' << this->latpar[1] << ',' <<
            this->latpar[2] << ',' << this->latpar[3] << ',' <<
            this->latpar[4] << ',' << this->latpar[5] << '\n';
        cout << "rmax=" << this->rmax << '\n';
    }
    // distreuse
    cout << "distreuse=" << distreuse << '\n';
    // tolcost
    cout << "tolcost=" << tolcost << '\n';
    // formula
    cout << "formula=" << this->formula.toString() << '\n';
    // radii
    if (!this->radii.empty())
    {
        cout << "radii=" << this->radii.toString() << '\n';
    }
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
    // promotefrac, promoterelax, demoterelax
    cout << "promotefrac=" << promotefrac << '\n';
    cout << "promoterelax=" << promoterelax << '\n';
    cout << "demoterelax=" << demoterelax << '\n';
    // ligasize, stopgame, seasontrials, trialsharing
    cout << "ligasize=" << ligasize << '\n';
    cout << "stopgame=" << stopgame << '\n';
    cout << "seasontrials=" << seasontrials << '\n';
    cout << "trialsharing=" << trialsharing << '\n';
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
    const char* parnames[] = {
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
        "crystal",
        "latpar",
        "rmax",
        "distreuse",
        "tolcost",
        "natoms",
        "formula",
        "radii",
        "fixed_atoms",
        "maxcputime",
        "rngseed",
        "promotefrac",
        "promoterelax",
        "demoterelax",
        "ligasize",
        "stopgame",
        "seasontrials",
	"trace",
        "trialsharing",
        "bangle_range",
        "max_dist",
        // obsolete ignored parameters
        "seed_clusters",
        "lookout_prob",
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
