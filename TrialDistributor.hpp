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

class TrialDistributor
{
    public:

	// class methods
	static TrialDistributor* create(const std::string& tp);
	static std::list<std::string> getTypes();
	static bool isType(const std::string& tp);

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
	const std::string& type()   { return _type; }
	virtual void share(int seasontrials) = 0;

    protected:

	// types
	typedef std::deque<double> BadnessHistory;

	// class data
	static const size_t histsize;

	// data members
	std::string _type;
	std::deque<BadnessHistory> lvbadlog;
	std::deque<double> fillrate;

	// protected methods

    private:

	// registry of derived trial distributors
	static std::map<std::string,TrialDistributor*> distributors;

	// methods
	virtual TrialDistributor* create() = 0;

};  // class TrialDistributor


class TrialDistributorEqual : public TrialDistributor
{
    public:

	// constructor and destructor
	TrialDistributorEqual()	    { _type = "equal"; }
	virtual ~TrialDistributorEqual() { }

	// methods
	virtual void share(int seasontrials);

    private:

	// methods
	virtual TrialDistributor* create()
	{
	    return new TrialDistributorEqual();
	}

};  // class TrialDistributorEqual


class TrialDistributorSize : public TrialDistributor
{
    public:

	// constructor and destructor
	TrialDistributorSize()	    { _type = "size"; }
	virtual ~TrialDistributorSize() { }

	// public methods
	virtual void share(int seasontrials);

    private:

	// methods
	virtual TrialDistributor* create()
	{
	    return new TrialDistributorSize();
	}

};  // class TrialDistributorSize


class TrialDistributorSuccess : public TrialDistributor
{
    public:

	// constructor and destructor
	TrialDistributorSuccess()   { _type = "success"; }
	virtual ~TrialDistributorSuccess() { }

	// public methods
	virtual void share(int seasontrials);

    private:

	// class data
	static const double success_weight;

	// methods
	virtual TrialDistributor* create()
	{
	    return new TrialDistributorSuccess();
	}

};  // class TrialDistributorSuccess


#endif	// TRIALDISTRIBUTOR_HPP_INCLUDED
