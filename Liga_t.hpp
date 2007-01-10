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

#include <set>
#include <vector>

#include "Division_t.hpp"
#include "RunPar_t.hpp"
class Molecule;

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
	inline bool stopFlag(int* flag)
	{
	    stopflag = flag;
	    return stopflag && *stopflag;
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

	// Private ethods
	void makeSeedClusters();
	PMOL updateWorldChamp();
	PMOL updateBestChamp();
	void printWorldChamp();
	void printLevelAverages();
	void saveOutStru();
	void saveFrames();
	void recordFramesTrace(const set<PMOL>& modified);
	void saveFramesTrace(const set<PMOL>& modified);
};

#endif		// LIGA_T_HPP_INCLUDED
