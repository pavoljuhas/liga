#ifndef BGAUTILS_HPP_INCLUDED
#define BGAUTILS_HPP_INCLUDED

#include <stdexcept>
#include <fstream>
#include <string>
#include <list>
#include <valarray>

template<class A, class E> inline void Assert(A assertion, E except)
{
    if (!assertion) throw except;
}

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

namespace BGA {
double CPUTime();
}   // namespace BGA

struct Counters_t
{
    long long penalty_calls;
    long long distance_calls;
    void PrintRunStats();
};
// global instance of Counters_t
namespace BGA {
extern Counters_t cnt;
}

double vdnorm(const std::valarray<double>&);
double vddot(const std::valarray<double>&, const std::valarray<double>&);
std::valarray<double> vdcross(const std::valarray<double>&,
			      const std::valarray<double>&);
std::valarray<double> vdrecipw0(const std::valarray<double>&);

#endif		// BGAUTILS_HPP_INCLUDED
