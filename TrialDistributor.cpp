/***********************************************************************
* Short Title: declaration of TrialDistributor and related classes
*
* Comments: 
*
* $Id$
* 
* <license text>
***********************************************************************/

#include "TrialDistributor.hpp"
#include "RunPar_t.hpp"

using namespace std;

RegisterSVNId TrialDistributor_cpp_id("$Id$");


////////////////////////////////////////////////////////////////////////
// class TrialDistributor
////////////////////////////////////////////////////////////////////////

// class data

map<std::string,TrialDistributor*> TrialDistributor::distributors;
const size_t TrialDistributor::histsize = 10;

// class methods

TrialDistributor* TrialDistributor::create(RunPar_t* rp)
{
    const string& tstp = rp->trials_sharing;
    if (!isType(tstp))
    {
	ostringstream emsg;
	emsg << "TrialDistributor '" << tstp << "' not defined.";
	throw runtime_error(emsg.str());
    }
    TrialDistributor* td = distributors[tstp]->create();
    // copy data from rp
    td->resize(rp->natoms + 1);
    td->tol_bad = rp->tol_bad;
    td->base_level = rp->base_level;
    return td;
}

list<string> TrialDistributor::getTypes()
{
    list<string> rv;
    typedef map<string,TrialDistributor*>::iterator Iter;
    for (Iter ii = distributors.begin(); ii != distributors.end(); ++ii)
    {
	rv.push_back(ii->first);
    }
    return rv;
}

bool TrialDistributor::isType(const string& tp)
{
    return distributors.count(tp);
}

// public methods

bool TrialDistributor::Register()
{
    if (distributors.count(_type))  return true;
    distributors[_type] = create();
    return true;
}

void TrialDistributor::setLevelBadness(size_t lv, double bd)
{
    BadnessHistory& hist = lvbadlog[lv];
    hist.push_back(bd);
    if (hist.size() > histsize)	    hist.pop_front();
}

void TrialDistributor::setLevelFillRate(size_t lv, double fr)
{
    fillrate[lv] = fr;
}

void TrialDistributor::resize(size_t sz)
{
    // clear history at every level
    lvbadlog.clear();
    lvbadlog.resize(sz);
    fillrate.resize(sz, 0.0);
    tshares.resize(sz, 0.0);
    top_level = sz - 1;
}

// protected methods

////////////////////////////////////////////////////////////////////////
// class TrialDistributorEqual
////////////////////////////////////////////////////////////////////////

void TrialDistributorEqual::share(int seasontrials)
{
    tshares = 0.0;
    if (base_level >= top_level)    return;
    tshares = 1.0 * seasontrials / (top_level - base_level);
    tshares[top_level] = 0.0;
}


////////////////////////////////////////////////////////////////////////
// class TrialDistributorSize
////////////////////////////////////////////////////////////////////////

void TrialDistributorSize::share(int seasontrials)
{
    tshares = 0.0;
    if (base_level >= top_level)    return;
    for (int lv = base_level; lv < top_level; ++lv)	tshares[lv] = lv;
    tshares[top_level] = 0.0;
    tshares *= seasontrials/tshares.sum();
}


////////////////////////////////////////////////////////////////////////
// class TrialDistributorSuccess
////////////////////////////////////////////////////////////////////////

// class data

void TrialDistributorSuccess::share(int seasontrials)
{
    tshares = 0.0;
    if (base_level >= top_level)    return;
    // calculate improvement ratio at each level
    valarray<double> scwt(0.0, size());
    for (int lv = base_level; lv <= top_level; ++lv)
    {
	BadnessHistory& hist = lvbadlog[lv];
	double tothwt = 0.0;
	int hwt = histsize - 1;
	for (int i = hist.size() - 1;  i > 0 && hwt > 0.0;  --i, --hwt)
	{
	    double improvement = (hist[i] < hist[i-1]) ?
		(hist[i-1] - hist[i])/tol_bad : 0;
	    scwt[lv] += hwt * improvement;
	    tothwt += hwt;
	}
	if (tothwt > 0.0)   scwt[lv] /= tothwt;
    }
    // split half of success weight to the levels below
    // for the last level move all success weitht to the former level
    for (int lo = 0, hi = 1; hi <= top_level; ++lo, ++hi)
    {
	double scshift = (hi < top_level) ? scwt[hi]/2 : scwt[hi];
	scwt[lo] += scshift;
	scwt[hi] -= scshift;
    }
    // calculate tshares
    // get size share weights
    valarray<double> szwt(0.0, size());
    for (int lv = base_level; lv < top_level; ++lv)	szwt[lv] = lv;
    szwt /= szwt.sum();
    // average success and size shares
    tshares = scwt + szwt;
    tshares[top_level] = 0.0;
    tshares *= seasontrials/tshares.sum();
}


////////////////////////////////////////////////////////////////////////
// registration of trial distributors
////////////////////////////////////////////////////////////////////////

namespace {

bool regTrialDistributorEqual = TrialDistributorEqual().Register();
bool regTrialDistributorSize = TrialDistributorSize().Register();
bool regTrialDistributorSuccess = TrialDistributorSuccess().Register();

}   // namespace


// End of TrialDistributor.cpp
