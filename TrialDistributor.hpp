/***********************************************************************
* Short Title: declaration of TrialDistributor and related classes
*
* Comments: 
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef TRIALDISTRIBUTOR_HPP_INCLUDED
#define TRIALDISTRIBUTOR_HPP_INCLUDED

#include <string>
#include <sstream>
#include <stdexcept>
#include <list>
#include <map>
#include <valarray>
#include <deque>
#include "RegisterSVNId.hpp"

namespace {
RegisterSVNId TrialDistributor_hpp_id = "$Id$";
}

class RunPar_t;

class TrialDistributor
{
    public:

	enum DistributorType { EQUAL, SIZE, SUCCESS };

	// class methods
	static TrialDistributor* create(RunPar_t* rp);
	static TrialDistributor* create(DistributorType tp);
	static std::list<std::string> getTypes();
	static bool isType(const std::string& tp);

	// constructor
	TrialDistributor()
	{
	    tol_bad = 1.0e-8;
	}
	// destructor
	virtual ~TrialDistributor() { }

	// data - calculated
	std::valarray<double> tshares;

	// function to register derived classes
	bool Register();

	// public methods
	void setLevelBadness(size_t lv, double bd);
	void setLevelFillRate(size_t lv, double fr);
	void resize(size_t sz);
	inline size_t size()	    { return lvbadlog.size(); }
	DistributorType type()	    { return _type; }
	const std::string& typeStr()	{ return _type_str; }
	virtual void share(int seasontrials) = 0;

    protected:

	// types
	typedef std::deque<double> BadnessHistory;

	// class data
	static const size_t histsize;

	// data members
	DistributorType _type;
	std::string _type_str;
	std::deque<BadnessHistory> lvbadlog;
	std::deque<double> fillrate;
	int base_level;
	int top_level;
	double tol_bad;

	// protected methods

    private:

	// registry of derived trial distributors
	static std::map<std::string,DistributorType> distributors;

};  // class TrialDistributor


class TrialDistributorEqual : public TrialDistributor
{
    public:

	// constructor and destructor
	TrialDistributorEqual()	: TrialDistributor()
	{
	    _type = EQUAL;
	    _type_str = "equal";
	}
	virtual ~TrialDistributorEqual() { }

	// methods
	virtual void share(int seasontrials);

};  // class TrialDistributorEqual


class TrialDistributorSize : public TrialDistributor
{
    public:

	// constructor and destructor
	TrialDistributorSize() : TrialDistributor()
	{
	    _type = SIZE;
	    _type_str = "size";
	}
	virtual ~TrialDistributorSize() { }

	// public methods
	virtual void share(int seasontrials);

};  // class TrialDistributorSize


class TrialDistributorSuccess : public TrialDistributor
{
    public:

	// constructor and destructor
	TrialDistributorSuccess() : TrialDistributor()
	{
	    _type = SUCCESS;
	    _type_str = "success";
	}
	virtual ~TrialDistributorSuccess() { }

	// public methods
	virtual void share(int seasontrials);

};  // class TrialDistributorSuccess


#endif	// TRIALDISTRIBUTOR_HPP_INCLUDED
