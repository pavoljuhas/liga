/***********************************************************************
* Short Title: figure optimum number of triangulations at each size
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include "TriangulationGuru.hpp"
#include "BGAutils.hpp"

using namespace std;

// Constructor

TriangulationGuru::TriangulationGuru(size_t max_atoms)
{
    check_size(max_atoms);
}


// Methods

vector<size_t> TriangulationGuru::shareTrials(size_t natoms, size_t n) const
{
    check_size(natoms);
    double p[3] = {
	equalbias + tginfos[natoms].tgcounts[0],
	equalbias + tginfos[natoms].tgcounts[1],
	equalbias + tginfos[natoms].tgcounts[2],
    };
    double ptot[3];
    partial_sum(p, p+3, ptot);
    std::vector<size_t> rv(3);
    switch (natoms)
    {
	case 1:
	    rv[0] = n;
	    break;
	case 2:
	    rv[0] = size_t( ceil(n*p[0]/ptot[1]) );
	    rv[1] = size_t( ceil(n*p[1]/ptot[1]) );
	    break;
	default:
	    for (size_t i = 0; i != 3; ++i)
	    {
		rv[i] = size_t( ceil(n*p[i]/ptot[2]) );
	    }
    }
    // make sure rv shares sum to n
    while (rv[0]+rv[1]+rv[2] > n)	*max(rv.begin(), rv.end()) -= 1;
    return rv;
}

size_t TriangulationGuru::trialsPerTime(size_t natoms, double cputime) const
{
    check_size(natoms);
    const double& rate = tginfos[natoms].trial_rate;
    double rv;
    rv = rate > 0.0 ? min(cputime/rate, double(maxtrials)) : mintrials;
    return size_t(rv);
}

void TriangulationGuru::tic()
{
    tictime = BGA::CPUTime();
}

void TriangulationGuru::toc(size_t natoms, size_t n)
{
    double toctime = BGA::CPUTime();
    check_size(natoms);
    tginfos[natoms].trial_rate = (toctime - tictime)/n;
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
