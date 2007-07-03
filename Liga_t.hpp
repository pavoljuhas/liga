/***********************************************************************
* Short Title: class Liga_t for competisions organization
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef LIGA_T_HPP_INCLUDED
#define LIGA_T_HPP_INCLUDED

#include <memory>
#include <set>
#include <vector>

#include "Division_t.hpp"
#include "RunPar_t.hpp"
#include "RegisterSVNId.hpp"

namespace {
RegisterSVNId Liga_t_hpp_id("$Id$");
}

class TrialDistributor;

class Liga_t : public std::vector<Division_t>
{
    public:

	// Data members
	int season;
	
	// Constructor and destructor
	Liga_t(RunPar_t* runpar);
	~Liga_t();

	// Public methods
	void prepare();
	void playSeason();
	void playLevel(size_t lo_level);
	inline bool stopFlag()
	{
	    return stopflag && *stopflag;
	}
	inline bool useStopFlag(int* flag)
	{
	    stopflag = flag;
	    return stopflag;
	}
	inline bool finished()
	{
	    isfinished = isfinished || stopflag && *stopflag ||
		solutionFound() || outOfTime();
	    return isfinished;
	}
	inline bool solutionFound()
	{
	    return world_champ && world_champ->Full() &&
		world_champ->NormBadness() < rp->tol_bad;
	}
	inline bool outOfTime()
	{
	    return rp->maxcputime > 0.0 && BGA::CPUTime() > rp->maxcputime;
	}
	void printFramesTrace();
	void printSummary();

    private:

	// Types
	typedef Molecule* PMOL;

	// Data members
	RunPar_t* rp;
	bool isfinished;
	int* stopflag;
	int base_level;
	PMOL world_champ, best_champ;
	std::auto_ptr<TrialDistributor> tdistributor;

	// Private methods
	int divSize(int level);
	void makeSeedClusters();
	void shareSeasonTrials();
	PMOL updateWorldChamp();
	PMOL updateBestChamp();
	void printWorldChamp();
	void printBestChamp();
	void printLevelAverages();
	void printTrialShares();
	void saveOutStru();
	void saveFrames();
	void recordFramesTrace(const set<PMOL>& modified);
	void saveFramesTrace(const set<PMOL>& modified);
};

#endif	// LIGA_T_HPP_INCLUDED
