/***********************************************************************
* Short Title: one division of the liga system
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef DIVISION_T_HPP_INCLUDED
#define DIVISION_T_HPP_INCLUDED

#include <vector>
#include "BGAlib.hpp"

class Division_t : public vector<Molecule*>
{
    private:

	// Types
	typedef Molecule* PMOL;
	// Data members
	size_t max_size;
	// Private methods
	static inline bool comp_PMOL_Badness(const PMOL lhs, const PMOL rhs)
	{
	    return lhs->Badness() < rhs->Badness();
	}

    public:

	// constructors
	Division_t(size_t s)  : vector<PMOL>(), max_size(s) { }
	Division_t(const Division_t& div0) :
	    vector<PMOL>(div0), max_size(div0.max_size) { }
	~Division_t()
	{
	    for (iterator ii = begin(); ii != end(); ++ii)
		delete *ii;
	}
	Division_t& operator= (const vector<PMOL>& div0)
	{
	    *this = div0;
	    return *this;
	}
	Division_t& operator= (const Division_t& div0)
	{
	    *this = vector<PMOL>(div0);
	    max_size = div0.max_size;
	    return *this;
	}
	int find_winner()
	{
	    // evaluate molecule fitness
	    valarray<double> vmfit(size());
	    double *pd = &vmfit[0];
	    for (iterator mi = begin(); mi != end(); ++mi, ++pd)
		*pd = (*mi)->NormBadness();
	    // then get the reciprocal value
	    vmfit = vdrecipw0(vmfit);
	    double *mfit = &vmfit[0];
	    int idx = random_wt_choose(1, mfit, size()).front();
	    return idx;
	}
	PMOL& winner()
	{
	    return at(find_winner());
	}
	int find_looser()
	{
	    // evaluate molecule fitness
	    valarray<double> vmbad(size());
	    double *pd = &vmbad[0];
	    for (iterator mi = begin(); mi != end(); ++mi, ++pd)
		*pd = (*mi)->NormBadness();
	    double *mbad = &vmbad[0];
	    int idx = random_wt_choose(1, mbad, size()).front();
	    return idx;
	}
	PMOL& looser()
	{
	    return at(find_looser());
	}
	PMOL& best()
	{
	    iterator pm = min_element(begin(), end(), comp_PMOL_Badness);
	    return *pm;
	}
	inline bool full() { return !(size() < max_size); }
	inline int fullsize() { return max_size; }
	inline double NormBadness()
	{
	    double total = 0.0;
	    for (iterator ii = begin(); ii != end(); ++ii)
	    {
		total += (*ii)->NormBadness();
	    }
	    return size() ? total / size() : 0.0;
	}

};

#endif		// DIVISION_T_HPP_INCLUDED
