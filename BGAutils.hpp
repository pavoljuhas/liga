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

template<class T> typename list<T>::iterator list_at(const list<T>& lst, int n);
double vdnorm(const valarray<double>&);
double vddot(const valarray<double>&, const valarray<double>&);
valarray<double> vdcross(const valarray<double>&, const valarray<double>&);
valarray<double> vdrecipw0(const valarray<double>&);

#endif		// BGAUTILS_HPP_INCLUDED
