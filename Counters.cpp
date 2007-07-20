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

#include "Counters.hpp"

using namespace std;
using namespace LIGA;

// Methods

double Counters::CPUTime()
{
    tms tbuf;
    times(&tbuf);
    return 1.0*tbuf.tms_utime/sysconf(_SC_CLK_TCK);
}

void Counters::printRunStats()
{
    char hostname[255];
    gethostname(hostname, 255);
    cout << "Run statistics:" << endl;
    cout << "dpenalty_calls = " << dpenalty_calls << endl;
    cout << "wpenalty_calls = " << wpenalty_calls << endl;
    cout << "R3_norm_calls = " << R3_norm_calls << endl;
    cout << "UserCPUtime = " << CPUTime() << 's' << endl;
    cout << "Host = " << hostname << endl;
}


// Global instance of Counters

Counters LIGA::counters;
