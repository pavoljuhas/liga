/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <sstream>
#include <stdexcept>
#include <boost/foreach.hpp>
#include "Crystal.hpp"
#include "Lattice.hpp"
#include "Atom_t.hpp"
#include "AtomCostCrystal.hpp"
#include "LigaUtils.hpp"
#include "AtomSequence.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// class Crystal
////////////////////////////////////////////////////////////////////////


// constructors

Crystal::Crystal() : Molecule()
{
    init();
}

Crystal::Crystal(const DistanceTable& dtab) : Molecule(dtab)
{
    init();
}

Crystal::Crystal(const Crystal& crs) : Molecule()
{
    init();
    *this = crs;
}

Crystal::~Crystal()
{ }

// public methods

Crystal& Crystal::operator=(const Crystal& crs)
{
    if (this == &crs)   return *this;
    // use base class copy
    Molecule& molthis = *this;
    molthis = crs;
    // copy Crystal specific data
    this->_lattice = crs._lattice;
    this->_rmin = crs._rmin;
    this->_rmax = crs._rmax;
    this->_cost = crs._cost;
    this->pmx_pair_counts = crs.pmx_pair_counts;
    this->_count_pairs = crs._count_pairs;
    this->_cost_data_cached = crs._cost_data_cached;
    return *this;
}


Molecule* Crystal::clone() const
{
    Molecule* pclone = new Crystal(*this);
    return pclone;
}


void Crystal::setDistanceTable(const DistanceTable& dtbl)
{
    DistanceTable dtblunique(dtbl.unique());
    this->Molecule::setDistanceTable(dtblunique);
    this->uncacheCostData();
}

void Crystal::setDistReuse(bool flag)
{
    if (!flag)
    {
        const char* emsg = "Crystal requires distreuse is true.";
        throw range_error(emsg);
    }
    this->Molecule::setDistReuse(flag);
}

void Crystal::setLattice(const Lattice& lat)
{
    _lattice.reset(new Lattice(lat));
    this->uncacheCostData();
}

const Lattice& Crystal::getLattice() const
{
    return *_lattice;
}

void Crystal::setRRange(double rmin, double rmax)
{
    _rmin = rmin;
    _rmax = rmax;
    this->uncacheCostData();
}

pair<double,double> Crystal::getRRange() const
{
    return make_pair(_rmin, _rmax);
}

// r-range extended by circum diameter of the primitive cell
pair<double,double> Crystal::getRExtent() const
{
    R3::Vector center(0.0, 0.0, 0.0);
    BOOST_FOREACH (Atom_t* ai, atoms)
    {
        center += ai->r;
    }
    center /= countAtoms();
    double max_offcenter = 0.0;
    BOOST_FOREACH (Atom_t* ai, atoms)
    {
        double ri = R3::norm(ai->r - center);
        if (ri > max_offcenter)     max_offcenter = ri;
    }
    // adjust for round-off errors, when not empty
    if (countAtoms() != 0)
    {
        const double epsd = sqrt(DOUBLE_EPS);
        max_offcenter = max_offcenter*(1.0 + epsd) + epsd;
    }
    double circum_diameter = 2*max_offcenter;
    return make_pair(_rmin - circum_diameter, _rmax + circum_diameter);
}

double Crystal::cost() const
{
    if (!this->_cost_data_cached)   recalculate();
    return this->Molecule::cost();
}


int Crystal::countPairs() const
{
    if (!this->_cost_data_cached)   recalculate();
    return this->_count_pairs;
}


void Crystal::recalculate() const
{
    // reset molecule
    this->ResetBadness();
    this->pmx_partial_costs.fill(0.0);
    this->pmx_pair_counts.fill(0);
    this->_count_pairs = 0;
    // reset all atoms
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
        seq.ptr()->ResetBadness();
    }
    AtomCostCrystal* atomcost;
    atomcost = static_cast<AtomCostCrystal*>(getAtomCostCalculator());
    // fill in diagonal elements
    R3::Vector zeros(0.0, 0.0, 0.0);
    pair<double,int> costcount;
    costcount = atomcost->pairCostCount(zeros, true);
    double diagpaircost = costcount.first;
    int diagpaircount = costcount.second;
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
        int idx = seq.ptr()->pmxidx;
        this->pmx_partial_costs(idx, idx) = diagpaircost;
        this->IncBadness(diagpaircost);
        seq.ptr()->IncBadness(diagpaircost);
        this->pmx_pair_counts(idx, idx) = diagpaircount;
        this->_count_pairs += diagpaircount;
    }
    // fill off-diagonal elements
    R3::Vector rc_dd;
    for (AtomSequence seq0(this); !seq0.finished(); seq0.next())
    {
        AtomSequence seq1 = seq0;
        if (!seq1.finished())   seq1.next();
        for (; !seq1.finished(); seq1.next())
        {
            int idx0 = seq0.ptr()->pmxidx;
            int idx1 = seq1.ptr()->pmxidx;
            rc_dd = seq1.ptr()->r - seq0.ptr()->r;
            costcount = atomcost->pairCostCount(rc_dd);
            // pmx is SymmetricMatrix, it is sufficient to make one assignment
            this->pmx_partial_costs(idx0, idx1) = costcount.first;
            this->IncBadness(costcount.first);
            double paircosthalf = costcount.first / 2.0;
            seq0.ptr()->IncBadness(paircosthalf);
            seq1.ptr()->IncBadness(paircosthalf);
            this->pmx_pair_counts(idx0, idx1) = costcount.second;
            this->_count_pairs += costcount.second;
        }
    }
    this->_cost_data_cached = true;
}


void Crystal::Clear()
{
    this->Molecule::Clear();
    uncacheCostData();
}

// protected methods

AtomCost* Crystal::getAtomCostCalculator() const
{
    static AtomCostCrystal the_acc(this);
    the_acc.resetFor(this);
    return &the_acc;
}

void Crystal::addNewAtomPairs(Atom_t* pa)
{
    if (!this->_cost_data_cached)   return;
    AtomCostCrystal* atomcost;
    atomcost = static_cast<AtomCostCrystal*>(getAtomCostCalculator());
    atomcost->eval(pa);
    // store partial costs
    const vector<double>& ptcs = atomcost->partialCosts();
    const vector<int>& pcnt = atomcost->pairCounts();
    assert(atoms.size() == ptcs.size());
    assert(atoms.size() == pcnt.size());
    for (AtomSequenceIndex seq(this); !seq.finished(); seq.next())
    {
	double paircost = ptcs[seq.idx()];
	int idx0 = pa->pmxidx;
	int idx1 = seq.ptr()->pmxidx;
	this->pmx_partial_costs(idx0, idx1) = paircost;
	double paircosthalf = paircost / 2.0;
	seq.ptr()->IncBadness(paircosthalf);
	pa->IncBadness(paircosthalf);
        int paircount = pcnt[seq.idx()];
        this->pmx_pair_counts(idx0, idx1) = paircount;
    }
    this->IncBadness(atomcost->totalCost());
    this->_count_pairs += atomcost->totalPairCount();
    // add self contribution:
    // calculates self cost and self pair count
    double diagpaircost;
    int diagpaircount;
    if (this->atoms.empty())
    {
        R3::Vector zeros(0.0, 0.0, 0.0);
        pair<double,int> costcount;
        costcount = atomcost->pairCostCount(zeros, true);
        diagpaircost = costcount.first;
        diagpaircount = costcount.second;
    }
    else
    {
        int idx1 = atoms.front()->pmxidx;
        diagpaircost = this->pmx_partial_costs(idx1, idx1);
        diagpaircount = this->pmx_pair_counts(idx1, idx1);
    }
    int idx0 = pa->pmxidx;
    this->pmx_partial_costs(idx0, idx0) = diagpaircost;
    pa->IncBadness(diagpaircost);
    this->IncBadness(diagpaircost);
    this->pmx_pair_counts(idx0, idx0) = diagpaircount;
    this->_count_pairs += diagpaircount;
    // take care of small round offs
    if (this->Badness() < NS_LIGA::eps_badness)	this->ResetBadness();
}


void Crystal::removeAtomPairs(Atom_t* pa)
{
    // remove associated pair costs and pair counts
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	int idx0 = pa->pmxidx;
	int idx1 = seq.ptr()->pmxidx;
	// remove pair costs
	double paircost = pmx_partial_costs(idx0, idx1);
	double paircosthalf = paircost/2.0;
	pa->DecBadness(paircosthalf);
	seq.ptr()->DecBadness(paircosthalf);
	this->DecBadness(paircost);
        // remove pair counts
        this->_count_pairs -= this->pmx_pair_counts(idx0, idx1);
    }
    if (this->Badness() < NS_LIGA::eps_badness) this->ResetBadness();
    assert(this->_count_pairs >= 0);
}


void Crystal::resizePairMatrices(int sz)
{
    int szcur = this->pmx_partial_costs.rows();
    if (sz < szcur)     return;
    Molecule::resizePairMatrices(sz);
    int sznew = this->pmx_partial_costs.rows();
    this->pmx_pair_counts.resize(sznew, sznew, 0);
}


// private methods

void Crystal::init()
{
    this->_rmin = 0.0;
    this->_rmax = 0.0;
    this->uncacheCostData();
    this->setDistReuse(true);
}


void Crystal::uncacheCostData()
{
    this->_cost_data_cached = false;
    if (countAtoms() == 0)
    {
        this->ResetBadness();
        this->_count_pairs = 0;
        this->pmx_pair_counts.fill(0);
        this->_cost_data_cached = true;
    }
}


// End of file
