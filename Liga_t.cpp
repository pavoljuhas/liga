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
#include "RunPar_t.hpp"
#include "TrialDistributor.hpp"

RegisterSVNId Liga_t_cpp_id("$Id$");


////////////////////////////////////////////////////////////////////////
// class Liga_t
////////////////////////////////////////////////////////////////////////

// Constructor and destructor

Liga_t::Liga_t(RunPar_t* runpar) :
    std::vector<Division_t>(), rp(runpar), stopflag(NULL)
{
    world_champ = NULL;
    best_champ = NULL;
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
    PMOL first_team = new Molecule(rp->mol);
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
//    const int* est_triang = lo_div->estimateTriangulations(tshares[lo_level]);
    const int* etg = lo_div->estimateTriangulations();
    advancing->Evolve(etg);
    lo_div->noteTriangulations(advancing);
    modified.insert(advancing);
    size_t hi_level = advancing->NAtoms();
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
    modified.insert(descending);
    // all set now so we can swap winner and looser
    hi_div->at(looser_idx) = advancing;
    lo_div->at(winner_idx) = descending;
    // make sure the original best cluster was much much better than
    // whatever left in the low division
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
    using namespace VerboseFlag;
    if (rp->verbose[AD])
    {
	cout << season;
	cout << " A " <<
	    lo_level << ' ' << adv_bad0 << ' ' <<
	    hi_level << ' ' << advancing->NormBadness() << "    ";
	cout << " D " <<
	    hi_level << ' ' << desc_bad0 << ' ' <<
	    lo_level << ' ' << descending->NormBadness() << '\n';
    }
    recordFramesTrace(modified);
    saveFramesTrace(modified);
    // update world champ so that the season can be cut short by stopFlag()
    if (advancing->Full())  updateWorldChamp();
}

void Liga_t::printFramesTrace()
{
    if (!rp->trace)	return;
    // needs clean up
    list<int>::iterator ii = best_champ->trace.begin();
    for (size_t n = 0; ii != best_champ->trace.end(); ++ii)
    {
	cout << ((n == 0) ? "TR " : " ") << *ii;
	if (++n == 3)
	{
	    cout << '\n';
	    n = 0;
	}
    }
    cout << endl;
}

void Liga_t::printSummary()
{

    if (solutionFound())    cout << "Solution found!!!\n\n";
    else if (outOfTime())   cout << "Exceeded maxcputime.\n\n";
    else if (stopFlag())    cout << "Simulation stopped, graceful death.\n\n";
    printFramesTrace();
    BGA::cnt.PrintRunStats();
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
	    while (mcore.NAtoms() > base_level) mcore.Pop(mcore.NAtoms() - 1);
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
    using namespace VerboseFlag;
    if (!rp->verbose[WC])   return;
    cout << season << " WC " << world_champ->NAtoms() << ' ' <<
	world_champ->NormBadness() << '\n';
}

void Liga_t::printBestChamp()
{
    using namespace VerboseFlag;
    if (!rp->verbose[BC])   return;
    cout << season << " BC " << best_champ->NAtoms() << ' '
	<< best_champ->NormBadness() << '\n';
}

void Liga_t::printLevelAverages()
{
    using namespace VerboseFlag;
    if (!rp->verbose[AV])   return;
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
    using namespace VerboseFlag;
    if (!rp->verbose[TS])    return;
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

    if (dontsave)   return;
    // start saving from the largest non-empty division
    for (reverse_iterator hi_div = rbegin(); hi_div != rend(); ++hi_div)
    {
	if (hi_div->empty())	continue;
	PMOL level_champ = hi_div->best();
	size_t level = level_champ->NAtoms();
	// save only if there is clear improvement
	if (!eps_lt(level_champ->NormBadness(), bestMNB[level]))    continue;
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
	// stop here if we are not saving all divisions
	if (!rp->saveall)   break;
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
    world_champ->WriteAtomEye(fname.c_str());
}

void Liga_t::recordFramesTrace(const set<PMOL>& modified)
{
    if (!rp->trace) return;
    Division_t::iterator mii;
    for (mii = at(base_level).begin(); mii != at(base_level).end(); ++mii)
    {
	(*mii)->trace.clear();
    }
    size_t level = base_level;
    for (iterator dv = begin() + base_level; dv != end(); ++dv, ++level)
    {
	size_t index = 0;
	for (mii = dv->begin(); mii != dv->end(); ++mii, ++index)
	{
	    PMOL mp = *mii;
	    if (!modified.count(mp) && !mp->trace.empty())  continue;
	    mp->trace.push_back(season);
	    mp->trace.push_back(level);
	    mp->trace.push_back(index);
	}
    }
}

void Liga_t::saveFramesTrace(const set<PMOL>& modified)
{
    static queue<TraceId_t> qtrace(rp->framestrace);
    while (!qtrace.empty() && season == qtrace.front().season)
    {
	TraceId_t tid = qtrace.front();
	if (tid.index >= int(at(tid.level).size()))	return;
	PMOL traced = at(tid.level).at(tid.index);
	if (!modified.count(traced))	return;
	qtrace.pop();
	int tno = rp->framestrace.size() - qtrace.size();
	ostringstream oss;
	oss << rp->frames << "." << tid.season << '.' << tno;
	string fname = oss.str();
	traced->WriteAtomEye(fname.c_str());
    }
}

// End of file
