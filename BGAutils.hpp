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

////////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////////

// constants

const double DOUBLE_MAX = std::numeric_limits<double>().max();
const double DOUBLE_EPS = std::numeric_limits<double>().epsilon();

namespace LIGA {
const double eps_badness = sqrt(DOUBLE_EPS);
};  // namespace LIGA

// functions

// round-off handling
bool eps_eq(const double& x, const double& y);
bool eps_gt(const double& x, const double& y);
bool eps_lt(const double& x, const double& y);

// valarray operations
double vdnorm(const std::valarray<double>&);
double vddot(const std::valarray<double>&,
        const std::valarray<double>&);
std::valarray<double> vdcross(const std::valarray<double>&,
			      const std::valarray<double>&);
// fitness
std::valarray<double> costToFitness(const std::valarray<double>&);
template <typename Container>
    Container costToFitness(const Container& vc);
template <typename Iterator>
    void costToFitnessInplace(Iterator first, Iterator last);

// file utilities
// similar to mkstemp(3)
std::ofstream& mktempofstream(std::ofstream& out, char *writefile);
bool read_header(std::istream& fid, std::string& header);
bool read_header(std::istream& fid);
template<typename T> bool read_data(std::istream& fid, std::vector<T>& v);


////////////////////////////////////////////////////////////////////////
// Definition for inline and template functions
////////////////////////////////////////////////////////////////////////

// round-off handling

inline bool eps_eq(const double& x, const double& y)
{
    return fabs(x-y) < LIGA::eps_badness;
}

inline bool eps_gt(const double& x, const double& y)
{
    return x > y + LIGA::eps_badness;
}

inline bool eps_lt(const double& x, const double& y)
{
    return x < y - LIGA::eps_badness;
}

// fitness

inline std::valarray<double> costToFitness(const std::valarray<double>& vc)
{
    std::valarray<double> vf(vc);
    double* first = &(vf[0]);
    double* last = &(vf[vf.size()]);
    costToFitnessInplace(first, last);
    return vf;
}

template <typename Container>
Container costToFitness(const Container& vc)
{
    Container vf(vc);
    costToFitnessInplace(vf.begin(), vf.end());
    return vf;
}

template <typename Iterator>
void costToFitnessInplace(Iterator first, Iterator last)
{
    for (Iterator ii = first; ii != last; ++ii)
    {
        *ii = 1.0 / (*ii + DOUBLE_EPS);
    }
}

// file utilities

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

#endif	// BGAUTILS_HPP_INCLUDED
