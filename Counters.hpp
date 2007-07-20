/***********************************************************************
* Short Title: LIGA statistics keeper
*
* Comments: class for keeping algorithm statistics
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef COUNTERS_HPP_INCLUDED
#define COUNTERS_HPP_INCLUDED

namespace LIGA {

class Counters
{
    public:

	// Data Members
	unsigned long long dpenalty_calls;
	unsigned long long wpenalty_calls;
	unsigned long long R3_norm_calls;

	// Methods
	double CPUTime();
	void printRunStats();

};

// Global instance of Counters
extern Counters counters;

}	// namespace LIGA

#endif	// COUNTERS_HPP_INCLUDED
