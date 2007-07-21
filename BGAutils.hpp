#ifndef BGAUTILS_HPP_INCLUDED
#define BGAUTILS_HPP_INCLUDED

#include <stdexcept>
#include <fstream>
#include <string>
#include <list>
#include <valarray>
#include <limits>
#include "RegisterSVNId.hpp"

namespace {
RegisterSVNId BGAutils_hpp_id("$Id$");
}

// constants
const double DOUBLE_MAX = std::numeric_limits<double>().max();

struct IOError : public std::runtime_error
{
    IOError (const std::string what_arg = "") :
	std::runtime_error(what_arg) { }
};

// similar to mkstemp(3)
std::ofstream& mktempofstream(std::ofstream& out, char *writefile);

template<typename T>
typename std::list<T>::iterator list_at(std::list<T>& lst, int n)
{
    typename std::list<T>::iterator ii;
    if (n <= lst.size()/2)
    {
	ii = lst.begin();
	advance(ii, n);
    }
    else
    {
	ii = lst.end();
	advance(ii, n-(int)lst.size());
    }
    return ii;
}

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

#endif		// BGAUTILS_HPP_INCLUDED
