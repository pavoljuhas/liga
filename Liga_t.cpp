/***********************************************************************
* Short Title: class Liga_t for competisions organization
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <queue>
#include "Liga_t.hpp"
#include "TraceId_t.hpp"
#include "Molecule.hpp"
#include "RunPar_t.hpp"
#include "TrialDistributor.hpp"

using namespace std;
using namespace NS_LIGA_VERBOSE_FLAG;

RegisterSVNId Liga_t_cpp_id("$Id$");


////////////////////////////////////////////////////////////////////////
// class Liga_t
////////////////////////////////////////////////////////////////////////

// class data

const vector<string> Liga_t::verbose_flags(verbose_flags_array,
        verbose_flags_array + VERBOSE_SIZE);

// class methods

bool Liga_t::isVerboseFlag(string flag)
{
    return count(verbose_flags.begin(), verbose_flags.end(), flag);
}

// Constructor and destructor

Liga_t::Liga_t(RunPar_t* runpar) :
    vector<Division_t>(), rp(runpar), stopflag(NULL)
{
    world_champ = NULL;
    best_champ = NULL;
    verbose = getDefaultVerbose();
    // this needs to be moved to RunPar_t later
    if (rp->args->ispar("verbose"))
    {
        setVerbose(ALL, false);
        vector<string>::iterator w;
        for (w = rp->verbose.begin(); w != rp->verbose.end(); ++w)
        {
            setVerbose(*w, true);
        }
    }
}

Liga_t::~Liga_t()
{
    delete best_champ;
    best_champ = NULL;
}

// Public methods

void Liga_t::prepare()
{
    season = 0;
    isfinished = false;
    clear();
    world_champ = NULL;
    delete best_champ;
    best_champ = NULL;
    tdistributor.reset( TrialDistributor::create(rp) );
    base_level = rp->base_level;
    // initialize divisions, primitive divisions have only 1 team
    Division_t::ndim = rp->ndim;
    for (int lev = 0; lev <= rp->natoms; ++lev)
    {
        int divsize = divSize(lev);
        push_back(Division_t(divsize, lev));
    }
    // put initial molecule to its division
    PMOL first_team = new Molecule(*rp->mol);
    cout << "Initial team" << endl;
    cout << season << " I " << first_team->NAtoms() << ' ' <<
	first_team->NormBadness() << '\n';
    at(first_team->NAtoms()).push_back(first_team);
    // fill lower divisions
    cout << "Filling lower divisions\n";
    for (int lev = first_team->NAtoms()-1; lev >= base_level; --lev)
    {
        PMOL parent_team = at(lev+1).back();
        PMOL lower_team = new Molecule(*parent_team);
        lower_team->Degenerate(1);
	cout << season << " L " << lower_team->NAtoms() << ' '
	    << lower_team->NormBadness() << endl;
        at(lev).push_back(lower_team);
    }
    makeSeedClusters();
    cout << "Done" << endl;
    updateWorldChamp();
    printWorldChamp();
    updateBestChamp();
    printBestChamp();
    cout << '\n' << "Starting the game ... now." << endl;
}

void Liga_t::playSeason()
{
    if (stopFlag())	return;
    ++season;
    shareSeasonTrials();
    for (size_t lo_level = base_level; lo_level < size() - 1; ++lo_level)
    {
	playLevel(lo_level);
	if (stopFlag())	    break;
    }
    printLevelAverages();
    updateWorldChamp();
    printWorldChamp();
    updateBestChamp();
    printBestChamp();
    cout.flush();
    saveOutStru();
    saveFrames();
}

void Liga_t::playLevel(size_t lo_level)
{
    iterator lo_div = begin() + lo_level;
    if (lo_div->empty())    return;
    // keep set of modified molecul for tracing
    set<PMOL> modified;
    // find winner
    int winner_idx = lo_div->find_winner();
    PMOL advancing = lo_div->at(winner_idx);
    if (advancing->NormBadness() >= rp->stopgame)    return;
    bool advancing_best = (advancing == lo_div->best());
    double adv_bad0 = advancing->NormBadness();
    if (!lo_div->full())
    {
	// save clone of advancing winner
	PMOL winner_clone = new Molecule(*advancing);
	lo_div->push_back(winner_clone);
    }
    // advance as far as possible
    const int* etg = lo_div->estimateTriangulations();
    advancing->Evolve(etg);
    lo_div->noteTriangulations(advancing);
    size_t hi_level = advancing->NAtoms();
    if (lo_level != hi_level)   modified.insert(advancing);
    iterator hi_div = begin() + hi_level;
    // fill intermediate empty divisions, loop will stop at non-empty lo_div
    for (iterator empty_div = hi_div; empty_div->empty(); --empty_div)
    {
	PMOL pioneer = new Molecule(*advancing);
	pioneer->Degenerate(hi_div - empty_div);
	empty_div->push_back(pioneer);
	modified.insert(pioneer);
    }
    // find looser
    int looser_idx = hi_div->find_looser();
    PMOL descending = hi_div->at(looser_idx);
    double desc_bad0 = descending->NormBadness();
    if (!hi_div->full())
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
    descending->Degenerate(hi_level - lo_level);
    if (lo_level != hi_level)   modified.insert(descending);
    // all set now so we can swap winner and looser
    hi_div->at(looser_idx) = advancing;
    lo_div->at(winner_idx) = descending;
    // make sure the original best cluster is preserved if it was much much
    // better than whatever left in the low division
    const double spoil_factor = 10.0;
    if ( advancing_best && eps_gt(lo_div->best()->NormBadness(),
		spoil_factor*adv_bad0) )
    {
	PMOL lo_looser = lo_div->at(lo_div->find_looser());
	*lo_looser = *advancing;
	for (size_t nlast = hi_level; nlast != lo_level;)
	    lo_looser->Pop(--nlast);
	modified.insert(lo_looser);
    }
    if (verbose[AD])
    {
	cout << season;
	cout << " A " <<
	    lo_level << ' ' << adv_bad0 << ' ' <<
	    hi_level << ' ' << advancing->NormBadness() << "    ";
	cout << " D " <<
	    hi_level << ' ' << desc_bad0 << ' ' <<
	    lo_level << ' ' << descending->NormBadness() << '\n';
    }
    recordFramesTrace(modified, lo_level);
    saveFramesTrace(modified, lo_level);
    // update world champ so that the season can be cut short by stopFlag()
    if (advancing->Full())  updateWorldChamp();
}

bool Liga_t::solutionFound() const
{
    return world_champ && world_champ->Full() &&
        world_champ->NormBadness() < rp->tol_bad;
}

void Liga_t::printFramesTrace() const
{
    if (!rp->trace)	return;
    // needs clean up
    cout << "Trace - season level natoms cost id:\n";
    list<TraceId_t>::iterator tii = best_champ->trace.begin();
    for (; tii != best_champ->trace.end(); ++tii)
    {
        cout << "TR " << tii->season <<
            ' ' << tii->level <<
            ' ' << tii->mol_natoms <<
            ' ' << tii->mol_norm_badness <<
            ' ' << tii->mol_id << '\n';
    }
    cout << endl;
}

void Liga_t::printSummary() const
{

    if (solutionFound())    cout << "Solution found!!!\n\n";
    else if (outOfTime())   cout << "Exceeded maxcputime.\n\n";
    else if (stopFlag())    cout << "Simulation stopped, graceful death.\n\n";
    printFramesTrace();
    Counter::printRunStats();
}

void Liga_t::setVerbose(VerboseFlag flag, bool value)
{
    if (flag < 0 || flag >= VERBOSE_SIZE)
    {
        ostringstream emsg;
        emsg << "Invalide verbose flag " << flag;
        throw invalid_argument(emsg.str());
    }
    verbose[flag] = value;
    if (flag == ALL)
    {
        fill(verbose.begin(), verbose.end(), value);
    }
}

void Liga_t::setVerbose(string flag, bool value)
{
    vector<string>::const_iterator vbsflag;
    vbsflag = find(verbose_flags.begin(), verbose_flags.end(), flag);
    if (vbsflag == verbose_flags.end())
    {
        ostringstream emsg;
        emsg << "Invalide verbose flag " << flag;
        throw invalid_argument(emsg.str());
    }
    int vbsindex = vbsflag - verbose_flags.begin();
    setVerbose(static_cast<VerboseFlag>(vbsindex), value);
}


const vector<bool>& Liga_t::getVerbose() const
{
    return verbose;
}

vector<bool> Liga_t::getDefaultVerbose()
{
    vector<bool> default_verbose(VERBOSE_SIZE, false);
    default_verbose[AD] = true;
    default_verbose[WC] = true;
    default_verbose[BC] = true;
    return default_verbose;
}

// Private methods

int Liga_t::divSize(int level)
{
    if (level < rp->base_level)	    return 0;
    if (level < 2)		    return 1;
    // default is rp->ligasize, but consult with seed_clusters
    int sz = rp->ligasize;
    for (   vector<SeedClusterInfo>::iterator scii = rp->seed_clusters.begin();
	    scii != rp->seed_clusters.end(); ++scii )
    {
	if (scii->level == level)   sz = scii->number;
    }
    return sz;
}

void Liga_t::makeSeedClusters()
{
    if (rp->seed_clusters.empty())  return;
    bool keep_evolve_jump = Molecule::evolve_jump;
    Molecule::evolve_jump = false;
    Molecule mcore( *(at(base_level).back()) );
    int keep_max_natoms = mcore.maxNAtoms();
    cout << "Generating seed clusters" << endl;
    for (vector<SeedClusterInfo>::iterator scii = rp->seed_clusters.begin();
	    scii != rp->seed_clusters.end(); ++scii )
    {
	mcore.setMaxNAtoms(scii->level);
	iterator seeded = begin() + scii->level;
	for (int nt = 0; nt < scii->trials; ++nt)
	{
	    // make sure mcore is at base_level
	    while (mcore.NAtoms() > base_level)
	    {
		mcore.Pop(mcore.NAtoms() - 1);
	    }
	    int addcnt = scii->level - base_level;
	    int ntrials = addcnt ? rp->seasontrials/addcnt + 1 : 0;
	    for (int k = 0; k < addcnt && !mcore.Full(); ++k)
	    {
		iterator lo_div = begin() + mcore.NAtoms();
		lo_div->assignTrials(ntrials);
		const int* etg = lo_div->estimateTriangulations();
		mcore.Evolve(etg);
	    }
	    if (!mcore.Full())  continue;
	    if (!seeded->full())
	    {
		seeded->push_back(new Molecule(mcore));
		continue;
	    }
	    // replace the worst cluster
	    PMOL replaced = seeded->at(seeded->find_looser());
	    *replaced = mcore;
	}
	// restore maxNAtoms
	for (size_t i = 0; i != seeded->size(); ++i)
	{
	    seeded->at(i)->setMaxNAtoms(keep_max_natoms);
	}
	PMOL seed_winner = seeded->at(seeded->find_winner());
	cout << season << " S " << seed_winner->NAtoms() << ' '
	    << seed_winner->NormBadness() << endl;
    }
    Molecule::evolve_jump = keep_evolve_jump;
}

void Liga_t::shareSeasonTrials()
{
    // copy level costs and fill rate to trials distributor
    for (size_t level = base_level; level != size(); ++level)
    {
	Division_t* lvdiv = &(at(level));
  	tdistributor->setLevelBadness(level, lvdiv->NormBadness());
	double fr = 1.0*lvdiv->size() / lvdiv->fullsize();
	tdistributor->setLevelFillRate(level, fr);
    }
    // share it
    tdistributor->share(rp->seasontrials);
    for (size_t level = base_level; level != size(); ++level)
    {
	Division_t* lvdiv = &(at(level));
	lvdiv->assignTrials(tdistributor->tshares[level]);
    }
    printTrialShares();
}

Molecule* Liga_t::updateWorldChamp()
{
    reverse_iterator ii;
    for (ii = rbegin(); ii != rend() && ii->empty(); ++ii)  { }
    world_champ = (ii != rend()) ? ii->best() : NULL;
    return world_champ;
}

Molecule* Liga_t::updateBestChamp()
{
    if (!world_champ)	return best_champ;
    if ( best_champ && world_champ->NAtoms() == best_champ->NAtoms() &&
	 !eps_lt(world_champ->NormBadness(), best_champ->NormBadness()) )
    {
	return best_champ;
    }
    // update needed here
    if (!best_champ)	best_champ = new Molecule(*world_champ);
    else		*best_champ = *world_champ;
    return best_champ;
}

void Liga_t::printWorldChamp()
{
    if (!verbose[WC])   return;
    cout << season << " WC " << world_champ->NAtoms() << ' ' <<
	world_champ->NormBadness() << '\n';
}

void Liga_t::printBestChamp()
{
    if (!verbose[BC])   return;
    cout << season << " BC " << best_champ->NAtoms() << ' '
	<< best_champ->NormBadness() << '\n';
}

void Liga_t::printLevelAverages()
{
    if (!verbose[AV])   return;
    size_t level = base_level;
    iterator lii = begin() + base_level;
    for (; lii != end() && !lii->empty(); ++lii, ++level)
    {
	cout << season << " AV " << level << ' ' <<
	    lii->NormBadness() << '\n';
    }
}

void Liga_t::printTrialShares()
{
    if (!verbose[TS])    return;
    for (size_t level = base_level; level < size(); ++level)
    {
	cout << season << " TS " << level <<
		' ' << tdistributor->tshares[level] << '\n';
    }
}

void Liga_t::saveOutStru()
{
    static int savecnt = 0;
    static valarray<double> bestMNB(DOUBLE_MAX, size());
    ++savecnt;
    bool dontsave = rp->outstru.empty() || this->empty() ||
	(!this->finished() && 0 < rp->saverate && savecnt < rp->saverate);
    // get out if there is nothing to do
    if (dontsave)   return;
    // find the largest non-empty division
    reverse_iterator save_div = rbegin();
    while (save_div != rend() && save_div->empty()) { ++save_div; }
    // save all lower division or stop after one save
    reverse_iterator stop_save;
    stop_save = (rp->saveall || save_div == rend()) ? rend() : save_div + 1;
    for (; save_div != stop_save; ++save_div)
    {
	// intermediate division may become empty
	if (save_div->empty())	continue;
	PMOL level_champ = save_div->best();
	size_t level = level_champ->NAtoms();
	// do not save empty structure
	if (level == 0)	break;
	// save only if there is clear improvement
	bool improved = eps_lt(level_champ->NormBadness(), bestMNB[level]);
	if (!improved)	continue;
	// something to save here
	savecnt = 0;
	bestMNB[level] = level_champ->NormBadness();
	// obtain file name
	string fname = rp->outstru;
	if (rp->saveall)
	{
	    ostringstream oss;
	    oss << rp->outstru << ".L" << level;
	    fname = oss.str();
	}
	level_champ->WriteFile(fname.c_str());
    }
}

void Liga_t::saveFrames()
{
    static struct {
	int cnt;
	PMOL champ;
	double level;
	double norm_badness;
    } saved = {0, NULL, 0, DOUBLE_MAX};
    bool dontsave = rp->frames.empty() || rp->framesrate == 0 ||
	++saved.cnt < rp->framesrate && !finished() || empty() ||
	world_champ == saved.champ && world_champ->NAtoms() == saved.level &&
	eps_eq(world_champ->NormBadness(), saved.norm_badness);
    if (dontsave)    return;
    // need to do something here
    saved.cnt = 0;
    saved.champ = world_champ;
    saved.level = world_champ->NAtoms();
    saved.norm_badness = world_champ->NormBadness();
    // build file name
    ostringstream oss;
    oss << rp->frames << '.' << season;
    string fname = oss.str();
    world_champ->WriteFile(fname.c_str());
}

void Liga_t::recordFramesTrace(set<PMOL>& modified, size_t lo_level)
{
    if (!rp->trace) return;
    Division_t::iterator mii;
    for (mii = at(base_level).begin(); mii != at(base_level).end(); ++mii)
    {
	(*mii)->trace.clear();
    }
    set<PMOL>::iterator modii;
    for (modii = modified.begin(); modii != modified.end(); ++modii)
    {
        PMOL mp = *modii;
        TraceId_t tid;
        tid.season = season;
        tid.level = lo_level;
        tid.mol_natoms = mp->NAtoms();
        tid.mol_norm_badness = mp->NormBadness();
        tid.mol_id = mp->id;
        mp->trace.push_back(tid);
    }
}

void Liga_t::saveFramesTrace(set<PMOL>& modified, size_t lo_level)
{
    bool dontsave = rp->frames.empty();
    if (dontsave)    return;
    static queue<TraceId_t> qtrace(rp->framestrace);
    while (!qtrace.empty())
    {
	TraceId_t tid = qtrace.front();
        if (season != tid.season)   break;
        if (long(lo_level) != tid.level)  break;
        PMOL traced = NULL;
        set<PMOL>::iterator modii = modified.begin();
        for (; !traced && modii != modified.end(); ++modii)
        {
            if (tid.mol_id == (*modii)->id)  traced = *modii;
        }
        if (!traced)    break;
        // traced molecule has been found here
	qtrace.pop();
	int tno = rp->framestrace.size() - qtrace.size();
	ostringstream oss;
	oss << rp->frames << "." << tid.season << '.' << tno;
	string fname = oss.str();
	traced->WriteFile(fname.c_str());
    }
}

// End of file
