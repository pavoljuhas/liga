#include <stdexcept>

#ifndef IOERROR_HPP_INCLUDED
#define IOERROR_HPP_INCLUDED

using namespace std;

struct IOError : public runtime_error
{
    IOError (const string what_arg = "") :
	runtime_error(what_arg) { }
};

#endif
