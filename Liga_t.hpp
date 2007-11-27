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
#include "Counter.hpp"
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
        static std::vector<bool> getDefaultVerbose();
        static void setVerboseVector(std::vector<bool>& vb,
                VerboseFlag flag, bool value=true);
        static void setVerboseVector(std::vector<bool>& vb,
                std::string flag, bool value=true);

        // instance data
	int season;

	// Constructor and destructor
	Liga_t(RunPar_t* runpar);
	~Liga_t();

	// Public methods
	void prepare();
	void playSeason();
	void playLevel(size_t lo_level);
	bool stopFlag() const;
	void useStopFlag(int* flag);
	bool finished() const;
	bool solutionFound() const;
	bool outOfTime() const;
	void printFramesTrace() const;
	void printSummary() const;
        void setVerbose(const std::vector<bool>& vb);
        const std::vector<bool>& getVerbose() const;

    private:

	// Types
	typedef Molecule* PMOL;

	// Data members
	RunPar_t* rp;
	int* stopflag;
	int base_level;
	PMOL world_champ, best_champ;
        bool printed_best_champ;
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
	void recordFramesTrace(std::set<PMOL>& modified, size_t lo_level);
	void saveFramesTrace(std::set<PMOL>& modified, size_t lo_level);
};

#endif	// LIGA_T_HPP_INCLUDED
