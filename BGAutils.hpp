#ifndef BGAUTILS_HPP_INCLUDED
#define BGAUTILS_HPP_INCLUDED

#include <stdexcept>
#include <list>
#include <valarray>

using namespace std;

struct IOError : public runtime_error
{
    IOError (const string what_arg = "") :
	runtime_error(what_arg) { }
};

template<class T> typename list<T>::iterator list_at(list<T>& lst, int n)
{
    typename list<T>::iterator ii;
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

struct Counters_t
{
    int penalty_calls;
    int distance_calls;
    void PrintRunStats();
};
// global instance of Counters_t
namespace BGA { extern Counters_t cnt; }

double vdnorm(const valarray<double>&);
double vddot(const valarray<double>&, const valarray<double>&);
valarray<double> vdcross(const valarray<double>&, const valarray<double>&);
valarray<double> vdrecipw0(const valarray<double>&);

#endif		// BGAUTILS_HPP_INCLUDED
