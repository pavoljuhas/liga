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

enum triangulation_type { LINEAR, PLANAR, SPATIAL };

class TriangulationGuru
{
    public:

	// Constructor
	TriangulationGuru(size_t max_atoms=0);

	// Methods
	std::vector<size_t> shareTrials(size_t natoms, size_t n) const;
	size_t trialsPerTime(size_t natoms, double cputime) const;
	// trial rate calibration
	void tic();
	void toc(size_t natoms, size_t n);
	// logging of successful triangulations
	void noteTriangulation(triangulation_type ttp, size_t natoms);
	void forgetTriangulations();

    private:

	// Constants
	static const double equalbias = 100.0;
	static const size_t mintrials = 128;
	static const size_t maxtrials = 65536;

	// Local helper struct
	struct TriangInfo
	{
	    size_t tgcounts[3];
	    double trial_rate;
	    TriangInfo() : trial_rate(0.0)
	    {
		tgcounts[0] = tgcounts[1] = tgcounts[2] = 0;
	    }
	};

	// Data members
	mutable std::vector<TriangInfo> tginfos;
	double tictime;

	// Methods
	void check_size(size_t natoms) const;

};

#endif		// TRIANGULATIONGURU_HPP_INCLUDED
