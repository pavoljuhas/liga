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
#include <string>

#include "Division_t.hpp"
#include "RunPar_t.hpp"
#include "RegisterSVNId.hpp"

namespace {
RegisterSVNId Liga_t_hpp_id("$Id$");
}

class TrialDistributor;

namespace NS_LIGA_VERBOSE_FLAG {

enum VerboseFlag { AD, WC, BC, AV, TS, ALL, VERBOSE_SIZE };

const std::string verbose_flags_array[VERBOSE_SIZE] = {
    "ad", "wc", "bc", "av", "ts", "all" };

}   // namespace NS_LIGA_VERBOSE_FLAG


class Liga_t : public std::vector<Division_t>
{
    public:

        // types
        typedef NS_LIGA_VERBOSE_FLAG::VerboseFlag VerboseFlag;

	// class data
        static const std::vector<std::string> verbose_flags;

        // class methods
        static bool isVerboseFlag(std::string flag);
        static std::vector<bool> getDefaultVerbose();

        // instance data
	int season;

	// Constructor and destructor
	Liga_t(RunPar_t* runpar);
	~Liga_t();

	// Public methods
	void prepare();
	void playSeason();
	void playLevel(size_t lo_level);
	inline bool stopFlag() const;
	inline bool useStopFlag(int* flag);
	inline bool finished();
	inline bool solutionFound() const;
	inline bool outOfTime() const;
	void printFramesTrace() const;
	void printSummary() const;
        void setVerbose(VerboseFlag flag, bool value=true);
        void setVerbose(std::string flag, bool value=true);
        const std::vector<bool>& getVerbose() const;

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
        std::vector<bool> verbose;

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

// class data

// class methods

// inline methods

inline bool Liga_t::stopFlag() const
{
    return stopflag && *stopflag;
}

inline bool Liga_t::useStopFlag(int* flag)
{
    stopflag = flag;
    return stopflag;
}

inline bool Liga_t::finished()
{
    isfinished = isfinished || stopflag && *stopflag ||
        solutionFound() || outOfTime();
    return isfinished;
}

inline bool Liga_t::solutionFound() const
{
    return world_champ && world_champ->Full() &&
        world_champ->NormBadness() < rp->tol_bad;
}

inline bool Liga_t::outOfTime() const
{
    return rp->maxcputime > 0.0 && BGA::CPUTime() > rp->maxcputime;
}

#endif	// LIGA_T_HPP_INCLUDED
