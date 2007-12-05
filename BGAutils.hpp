#ifndef BGAUTILS_HPP_INCLUDED
#define BGAUTILS_HPP_INCLUDED

#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
#include <limits>
#include "RegisterSVNId.hpp"

namespace {
RegisterSVNId BGAutils_hpp_id("$Id$");
}

// constants
const double DOUBLE_MAX = std::numeric_limits<double>().max();
const double DOUBLE_EPS = std::numeric_limits<double>().epsilon();

namespace LIGA {

const double eps_badness = sqrt(DOUBLE_EPS);

};

inline bool eps_eq(double x, double y)
{
    return fabs(x-y) < LIGA::eps_badness;
}

inline bool eps_gt(double x, double y)
{
    return x > y + LIGA::eps_badness;
}

inline bool eps_lt(double x, double y)
{
    return x < y - LIGA::eps_badness;
}

// similar to mkstemp(3)
std::ofstream& mktempofstream(std::ofstream& out, char *writefile);

double vdnorm(const std::valarray<double>&);
double vddot(const std::valarray<double>&, const std::valarray<double>&);
std::valarray<double> vdcross(const std::valarray<double>&,
			      const std::valarray<double>&);

// template function for zero-safe calculation of reciprocal vector
const double zero_reciprocal_gain = 10.0;
template <typename T>
T recipw0(const T& v, double zerogain=zero_reciprocal_gain)
{
    // create container with return values
    T rv = v;
    double* rvfirst = &rv[0];
    double* rvlast = &rv[rv.size()];
    double min_positive = DOUBLE_MAX;
    for (double* p = rvfirst; p != rvlast; ++p)
    {
	double& rvi = *p;
	if (0.0 < rvi && rvi < min_positive)	min_positive = rvi;
    }
    if (min_positive == DOUBLE_MAX)	min_positive = 1.0;
    double reczero = zerogain*1.0/min_positive;
    // calculate reciprocal values
    for (double* p = rvfirst; p != rvlast; ++p)
    {
	double& rvi = *p;
	rvi = (rvi != 0.0) ? 1.0/rvi : reczero;
    }
    return rv;
}

bool read_header(std::istream& fid, std::string& header);
bool read_header(std::istream& fid);
template<typename T> bool read_data(std::istream& fid, std::vector<T>& v);

////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////

// template functions
template<typename T>
bool read_data(std::istream& fid, std::vector<T>& v)
{
    // prepare v
    T x;
    while (fid >> x)
    {
	v.push_back(x);
    }
    return !(fid.rdstate() & std::ios::badbit);
}



#endif		// BGAUTILS_HPP_INCLUDED
