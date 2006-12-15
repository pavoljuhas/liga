/***********************************************************************
* Short Title: figure optimum number of triangulations at each size
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <stdexcept>
#include <numeric>
#include <iostream>

#include "BGAutils.hpp"
#include "TriangulationGuru.hpp"

using namespace std;

// Global data members
size_t TriangulationGuru::ndim = 3;

// Constructor

TriangulationGuru::TriangulationGuru(size_t max_atoms)
{
    fill(est_trials, est_trials+3, size_t(0));
    check_size(max_atoms);
}


// Methods

void TriangulationGuru::shareTrials(size_t natoms, size_t n) const
{
    check_size(natoms);
    double p[3] = {
	10 + tginfos[natoms].tgcounts[0],
	50 + tginfos[natoms].tgcounts[1],
	1000 + tginfos[natoms].tgcounts[2],
    };
    double ptot[3];
    partial_sum(p, p+3, ptot);
    fill(est_trials, est_trials+3, 0);
    // how many dimensions should be searched?
    size_t nd = min(ndim, natoms);
    for (size_t i = 0; i != nd; ++i)
    {
	est_trials[i] = size_t( ceil(n*p[i]/ptot[nd-1]) );
    }
    // make sure est_trials shares sum to n
    for (size_t etsum = accumulate(est_trials, est_trials + 3, 0); etsum > n; )
    {
	*max_element(est_trials, est_trials + 3) -= 1;
	etsum -= 1;
    }
cout << natoms << " n = " << n << "  est_trials = { ";
copy(est_trials, est_trials + 3, ostream_iterator<size_t,char>(cout," "));
cout << "}\n";
}

size_t TriangulationGuru::trialsPerDistanceCalls(
	size_t natoms, long long dcalls) const
{
    check_size(natoms);
    double rate = tginfos[natoms].trial_rate();
    double rv;
    rv = min(dcalls*rate, double(maxtrials));
    rv = max(rv, double(mintrials));
    return size_t(rv);
}

void TriangulationGuru::tic()
{
    tic_dcalls = BGA::cnt.distance_calls;
}

void TriangulationGuru::toc(size_t natoms, size_t n)
{
    if (tic_dcalls < 0)	    return;
    check_size(natoms);
    tginfos[natoms].tot_trials += n;
    tginfos[natoms].tot_dcalls += BGA::cnt.distance_calls - tic_dcalls;
    tic_dcalls = -1;
}

void TriangulationGuru::noteTriangulation(triangulation_type ttp,
	size_t natoms)
{
    check_size(natoms);
    tginfos[natoms].tgcounts[ttp] += 1;
}

void TriangulationGuru::forgetTriangulations()
{
    tginfos.clear();
}

void TriangulationGuru::check_size(size_t natoms) const
{
    if (tginfos.size() >= natoms + 1)	    return;
    tginfos.resize(natoms + 1);
}

// End of file
