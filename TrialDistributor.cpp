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

using namespace std;

RegisterSVNId TrialDistributor_cpp_id("$Id$");


////////////////////////////////////////////////////////////////////////
// class TrialDistributor
////////////////////////////////////////////////////////////////////////

// class data

map<std::string,TrialDistributor*> TrialDistributor::distributors;
const size_t TrialDistributor::histsize = 10;

// class methods

TrialDistributor* TrialDistributor::create(const string& tp)
{
    if (!isType(tp))
    {
	ostringstream emsg;
	emsg << "TrialDistributor '" << tp << "' not defined.";
	throw runtime_error(emsg.str());
    }
    return distributors[tp]->create();
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
}

// protected methods

////////////////////////////////////////////////////////////////////////
// class TrialDistributorEqual
////////////////////////////////////////////////////////////////////////

void TrialDistributorEqual::share(int seasontrials)
{
    tshares = 0.0;
    if (size() < 2)	return;
    tshares = 1.0 * seasontrials / (size() - 1);
    tshares[size() - 1] = 0.0;
}


////////////////////////////////////////////////////////////////////////
// class TrialDistributorSize
////////////////////////////////////////////////////////////////////////

void TrialDistributorSize::share(int seasontrials)
{
    tshares = 0.0;
    if (size() < 2)	return;
    for (size_t lv = 0; lv < size() - 1; ++lv)	tshares[lv] = lv;
    tshares *= seasontrials/tshares.sum();
    tshares[size() - 1] = 0.0;
}


////////////////////////////////////////////////////////////////////////
// class TrialDistributorSuccess
////////////////////////////////////////////////////////////////////////

// class data

const double TrialDistributorSuccess::tol_bad_scale = 0.01;

void TrialDistributorSuccess::share(int seasontrials)
{
    tshares = 0.0;
    if (size() < 2)	return;
    const double eps_tol_bad = tol_bad_scale * tol_bad;
    // calculate improvement ratio at each level
    valarray<double> scwt(0.0, size());
    for (size_t lv = 0; lv < size(); ++lv)
    {
	BadnessHistory& hist = lvbadlog[lv];
	double tothwt = 0.0;
	int hwt = histsize - 1;
	for (int i = hist.size() - 1;  i > 0 && hwt > 0.0;  --i, --hwt)
	{
	    double improvement = hist[i] + eps_tol_bad < hist[i-1] ?
		(hist[i-1] - hist[i])/(hist[i-1] + eps_tol_bad) : 0;
	    scwt[lv] += hwt * improvement;
	    tothwt += hwt;
	}
	if (tothwt > 0.0)   scwt[lv] /= tothwt;
    }
    // split half of success weight to the levels below
    // for the last level move all success weitht to the former level
    for (size_t lo = 0, hi = 1; hi < size(); ++lo, ++hi)
    {
	double scshift = (hi + 1 < size()) ? scwt[hi]/2 : scwt[hi];
	scwt[lo] += scshift;
	scwt[hi] -= scshift;
    }
    // calculate tshares
    // average with equal shares to work around zero scwt
    tshares = scwt + 1.0/(size() - 1);
    tshares[size() - 1] = 0.0;
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
