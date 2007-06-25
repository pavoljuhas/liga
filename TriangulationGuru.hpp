/***********************************************************************
* Short Title: figure optimum number of triangulations at each size
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef TRIANGULATIONGURU_HPP_INCLUDED
#define TRIANGULATIONGURU_HPP_INCLUDED

#include <vector>

class TriangulationGuru
{
    public:

	// Data members
	// Calculated:
	mutable size_t est_trials[3];
	// Constructor
	TriangulationGuru(size_t max_atoms=0);

	// Methods
	void shareTrials(size_t natoms, size_t n) const;
	size_t trialsPerDistanceCalls(size_t natoms, long long dcalls) const;
	// trial rate calibration
	void tic();
	void toc(size_t natoms, size_t n);
	// logging of successful triangulations
	void noteTriangulation(triangulation_type ttp, size_t natoms);
	void forgetTriangulations();

	// Constants
	static const size_t mintrials = 64;
	static const size_t maxtrials = 16384;
	static size_t ndim;

    private:

	// Local helper struct
	struct TriangInfo
	{
	    size_t tgcounts[3];
	    long long tot_trials;
	    long long tot_dcalls;
	    inline double trial_rate()
	    {
		return tot_dcalls ? double(tot_trials)/tot_dcalls : 0.0;
	    }
	    TriangInfo() : tot_trials(0), tot_dcalls(0)
	    {
		std::fill(tgcounts, tgcounts+3, 0);
	    }
	};

	// Data members
	mutable std::vector<TriangInfo> tginfos;
	long long tic_dcalls;

	// Methods
	void check_size(size_t natoms) const;

};

#endif	// TRIANGULATIONGURU_HPP_INCLUDED
