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
#include "LigaUtils.hpp"
#include "TraceId_t.hpp"
#include "Molecule.hpp"
#include "RunPar_t.hpp"
#include "TrialDistributor.hpp"

using namespace std;
using namespace NS_LIGA_VERBOSE_FLAG;

////////////////////////////////////////////////////////////////////////
// class Liga_t
////////////////////////////////////////////////////////////////////////

// class data

const vector<string> Liga_t::verbose_flags(verbose_flags_array,
        verbose_flags_array + VERBOSE_SIZE);

// class methods

vector<bool> Liga_t::getDefaultVerbose()
{
    vector<bool> default_verbose(VERBOSE_SIZE, false);
    default_verbose[WC] = true;
    default_verbose[BC] = true;
    return default_verbose;
}

void Liga_t::setVerboseVector(vector<bool>& vb, VerboseFlag flag, bool value)
{
    assert(vb.size() == VERBOSE_SIZE);
    if (flag < 0 || flag >= VERBOSE_SIZE)
    {
        ostringstream emsg;
        emsg << "Invalid verbose flag '" << flag << "'.";
        throw invalid_argument(emsg.str());
    }
    if (flag == ALL)    vb.assign(VERBOSE_SIZE, value);
    else                vb[flag] = value;
}

void Liga_t::setVerboseVector(vector<bool>& vb, string flag, bool value)
{
    assert(vb.size() == VERBOSE_SIZE);
    vector<string>::const_iterator vbsflag;
    vbsflag = find(verbose_flags.begin(), verbose_flags.end(), flag);
    if (vbsflag == verbose_flags.end())
    {
        ostringstream emsg;
        emsg << "Invalid verbose flag " << flag;
        throw invalid_argument(emsg.str());
    }
    int vbsindex = vbsflag - verbose_flags.begin();
    VerboseFlag vbf = static_cast<VerboseFlag>(vbsindex);
    setVerboseVector(vb, vbf, value);
}

// Constructor and destructor

Liga_t::Liga_t(RunPar_t* runpar) :
    vector<Division_t>(), rp(runpar), stopflag(NULL)
{
    world_champ = NULL;
    setVerbose(rp->verbose);
}

Liga_t::~Liga_t()
{
}

// Public methods

void Liga_t::prepare()
{
    season = 0;
    clear();
    this->world_champ = NULL;
    this->best_champ.reset(NULL);
    this->printed_best_champ = false;
    this->tdistributor.reset( TrialDistributor::create(rp) );
    this->base_level = rp->base_level;
    // initialize divisions, primitive divisions have only 1 team
    Division_t::ndim = rp->ndim;
    for (int lev = 0; lev <= rp->mol->getMaxAtomCount(); ++lev)
    {
        int divsize = divSize(lev);
        push_back(Division_t(divsize, lev));
    }
    // put initial molecule to its division
    PMOL first_team = rp->mol->clone();
    cout << "Initial team" << endl;
    cout << season << " I " << first_team->countAtoms() << ' ' <<
	first_team->cost() << '\n';
    at(first_team->countAtoms()).push_back(first_team);
    // fill lower divisions
    cout << "Filling lower divisions\n";
    for (int lev = first_team->countAtoms()-1; lev >= base_level; --lev)
    {
        PMOL parent_team = at(lev+1).back();
        PMOL lower_team = parent_team->clone();
        lower_team->Degenerate(1);
	cout << season << " L " << lower_team->countAtoms() << ' '
	    << lower_team->cost() << endl;
        at(lev).push_back(lower_team);
    }
    cout << "Done" << endl;
    updateWorldChamp();
    printWorldChamp();
    updateBestChamp();
    printBestChamp();
    saveOutStru();
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
    if (advancing->cost() >= rp->stopgame)    return;
    bool advancing_best = (advancing == lo_div->best());
    double adv_bad0 = advancing->cost();
    if (!lo_div->full())
    {
	// save clone of advancing winner
	PMOL winner_clone = advancing->clone();
	lo_div->push_back(winner_clone);
    }
    // advance as far as possible
    const int* etg = lo_div->estimateTriangulations();
    const pair<int*,int*> acc_tot = advancing->Evolve(etg);
    lo_div->noteTriangulations(acc_tot);
    size_t hi_level = advancing->countAtoms();
    if (lo_level != hi_level)   modified.insert(advancing);
    iterator hi_div = begin() + hi_level;
    // fill intermediate empty divisions, loop will stop at non-empty lo_div
    for (iterator empty_div = hi_div; empty_div->empty(); --empty_div)
    {
	PMOL pioneer = advancing->clone();
	pioneer->Degenerate(hi_div - empty_div);
	empty_div->push_back(pioneer);
	modified.insert(pioneer);
    }
    // find looser
    int looser_idx = hi_div->find_looser();
    PMOL descending = hi_div->at(looser_idx);
    double desc_bad0 = descending->cost();
    if (!hi_div->full())
    {
	// save clone of descending looser
	PMOL looser_clone = descending->clone();
	hi_div->push_back(looser_clone);
    }
    // clone winner if he made a good advance
    if (eps_gt(descending->cost(), advancing->cost()))
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
    if ( advancing_best && eps_gt(lo_div->best()->cost(),
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
	    hi_level << ' ' << advancing->cost() << "    ";
	cout << " D " <<
	    hi_level << ' ' << desc_bad0 << ' ' <<
	    lo_level << ' ' << descending->cost() << '\n';
    }
    recordFramesTrace(modified, lo_level);
    saveFramesTrace(modified, lo_level);
    // update world champ so that the season can be cut short by stopFlag()
    if (advancing->full())  updateWorldChamp();
}

bool Liga_t::stopFlag() const
{
    return stopflag && *stopflag;
}

void Liga_t::useStopFlag(int* flag)
{
    stopflag = flag;
}

bool Liga_t::finished() const
{
    bool isfinished = stopFlag() || solutionFound() || outOfTime();
    return isfinished;
}

bool Liga_t::solutionFound() const
{
    return world_champ && world_champ->full() &&
        world_champ->cost() < rp->tolcost;
}

bool Liga_t::outOfTime() const
{
    return rp->maxcputime > 0.0 && Counter::CPUTime() > rp->maxcputime;
}

void Liga_t::printFramesTrace() const
{
    if (!rp->trace)	return;
    // needs clean up
    cout << "Trace - season level natoms cost id:\n";
    list<TraceId_t>::iterator tii = this->best_champ->trace.begin();
    for (; tii != this->best_champ->trace.end(); ++tii)
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

void Liga_t::setVerbose(const vector<bool>& vb)
{
    assert(vb.size() == VERBOSE_SIZE);
    verbose = vb;
}

const vector<bool>& Liga_t::getVerbose() const
{
    return verbose;
}

// Private methods

int Liga_t::divSize(int level)
{
    // start with the default value
    int sz = rp->ligasize;
    // ignore divisions below base_level
    if (level < rp->base_level)	        sz = 0;
    // it is enough to have one team at base_level
    else if (level == rp->base_level)	sz = 1;
    // and also at levels 0 and 1
    else if (level < 2)		        sz = 1;
    return sz;
}


void Liga_t::shareSeasonTrials()
{
    // copy level costs and fill rate to trials distributor
    for (size_t level = base_level; level != size(); ++level)
    {
	Division_t* lvdiv = &(at(level));
  	tdistributor->setLevelBadness(level, lvdiv->averageCost());
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

void Liga_t::updateBestChamp()
{
    bool hasnewchamp =
        (this->world_champ && !this->best_champ.get()) ||
        this->world_champ->countAtoms() > this->best_champ->countAtoms() ||
        eps_lt(this->world_champ->cost(), this->best_champ->cost());
    if (hasnewchamp)
    {
        this->injectOverlapMinimization();
        this->best_champ.reset(this->world_champ->clone());
        this->printed_best_champ = false;
    }
}

void Liga_t::printWorldChamp()
{
    if (!verbose[WC])   return;
    cout << season << " WC " << world_champ->countAtoms() << ' ' <<
	world_champ->cost() << '\n';
}

void Liga_t::printBestChamp()
{
    bool dontprint = !verbose[BC] || (printed_best_champ && !finished());
    if (dontprint)	return;
    cout << this->season << " BC " << 
        this->best_champ->countAtoms() << ' '
	<< this->best_champ->cost() << endl;
    this->printed_best_champ = true;
}

void Liga_t::printLevelAverages()
{
    if (!verbose[AV])   return;
    size_t level = base_level;
    iterator lii = begin() + base_level;
    for (; lii != end() && !lii->empty(); ++lii, ++level)
    {
	cout << season << " AV " << level << ' ' <<
	    lii->averageCost() << '\n';
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
    static valarray<double> bestMcost(DOUBLE_MAX, size());
    ++savecnt;
    bool dontsave = rp->outstru.empty() || this->empty() ||
	(!this->finished() && (rp->saverate == 0 || savecnt < rp->saverate));
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
	size_t level = level_champ->countAtoms();
	// do not save empty structure
	if (level == 0)	break;
	// save only if there is clear improvement
	bool improved = eps_lt(level_champ->cost(), bestMcost[level]);
	if (!improved)	continue;
	// something to save here
	savecnt = 0;
	bestMcost[level] = level_champ->cost();
	// construct file name
	ostringstream fname;
        fname << rp->outstru;
	if (rp->saveall)
        {
	    fname << ".L" << level;
	}
	level_champ->WriteFile(fname.str().c_str());
    }
}

void Liga_t::saveFrames()
{
    static struct {
	int cnt;
	PMOL champ;
	double level;
	double cost;
    } saved = {0, NULL, 0, DOUBLE_MAX};
    bool dontsave = rp->frames.empty() || rp->framesrate == 0 ||
	(++saved.cnt < rp->framesrate && !finished()) || empty() ||
        (world_champ == saved.champ &&
         world_champ->countAtoms() == saved.level &&
         eps_eq(world_champ->cost(), saved.cost));
    if (dontsave)    return;
    // need to do something here
    saved.cnt = 0;
    saved.champ = world_champ;
    saved.level = world_champ->countAtoms();
    saved.cost = world_champ->cost();
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
        tid.mol_natoms = mp->countAtoms();
        tid.mol_norm_badness = mp->cost();
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

void Liga_t::injectOverlapMinimization()
{
    world_champ->DownhillOverlapMinimization();
    auto_ptr<Molecule> injector(world_champ->clone());
    while (injector->countAtoms() > base_level + 1)
    {
        injector->Degenerate(1);
        int level = injector->countAtoms();
        assert(!at(level).empty());
        int tgt_idx = at(level).find_looser();
	PMOL tgt_mol = at(level).at(tgt_idx);
        *tgt_mol = *injector;
    }
}

// End of file
