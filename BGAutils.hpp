#include <stdexcept>

#ifndef BGAUTILS_HPP_INCLUDED
#define BGAUTILS_HPP_INCLUDED

using namespace std;

struct IOError : public runtime_error
{
    IOError (const string what_arg = "") :
	runtime_error(what_arg) { }
};

#endif
