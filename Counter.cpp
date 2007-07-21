/***********************************************************************
* Short Title: LIGA statistics keeper
*
* Comments: class for keeping algorithm statistics
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <iostream>
#include <unistd.h>
#include <sys/times.h>

#include "Counter.hpp"

using namespace std;

// class data - private

Counter::CounterStorage Counter::storage;


// class methods

Counter* Counter::getCounter(string name)
{
    if (!storage.count(name))
    {
        storage[name] = new Counter(name);
    }
    return storage[name];
}

double Counter::CPUTime()
{
    tms tbuf;
    times(&tbuf);
    return 1.0*tbuf.tms_utime/sysconf(_SC_CLK_TCK);
}

void Counter::printRunStats()
{
    char hostname[255];
    gethostname(hostname, 255);
    cout << "Run statistics:\n";
    CounterStorage::iterator ii;
    for (ii = storage.begin(); ii != storage.end(); ++ii)
    {
        Counter& cii = *(ii->second);
        cout << cii << '\n';
    }
    cout << "UserCPUtime = " << CPUTime() << "s\n";
    cout << "Host = " << hostname << endl;
}


// constructor - private

Counter::Counter(string name) : _name(name)
{
    reset();
}


// non-member operators

ostream& operator<<(ostream& os, const Counter& cnt)
{
    os << cnt.name() << " = " << cnt.value();
    return os;
}

////////////////////////////////////////////////////////////////////////
// definitions for Counter::CounterStorage
////////////////////////////////////////////////////////////////////////

Counter::CounterStorage::~CounterStorage()
{
    for (iterator ii = begin(); ii != end(); ++ii)
    {
        delete ii->second;
        ii->second = NULL;
    }
}

// End of file
