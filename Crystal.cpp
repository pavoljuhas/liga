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
using namespace NS_LIGA;

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


Crystal& Crystal::operator=(const Molecule& mol)
{
    const Crystal& crs = dynamic_cast<const Crystal&>(mol);
    return operator=(crs);
}


Crystal& Crystal::operator=(const Crystal& crs)
{
    if (this == &crs)   return *this;
    // execute base assignment
    Molecule::operator=(crs);
    // copy Crystal specific data
    this->_lattice = crs._lattice;
    this->_rmin = crs._rmin;
    this->_rmax = crs._rmax;
    this->_cost = crs._cost;
    this->_full_distance_table = crs._full_distance_table;
    this->pmx_pair_counts = crs.pmx_pair_counts;
    this->_count_pairs = crs._count_pairs;
    this->_cost_data_cached = crs._cost_data_cached;
    this->_rextent = crs._rextent;
    this->_rextent_cached = crs._rextent_cached;
    return *this;
}


Molecule* Crystal::clone() const
{
    Molecule* pclone = new Crystal(*this);
    return pclone;
}


void Crystal::setDistanceTable(const DistanceTable& dtbl)
{
    vector<double> dtbl_unique = dtbl.unique();
    this->_full_distance_table.reset(new DistanceTable(dtbl_unique));
    this->cropDistanceTable();
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
    this->uncacheRExtent();
}

const Lattice& Crystal::getLattice() const
{
    return *_lattice;
}

void Crystal::setRmax(double rmax)
{
    this->_rmax = rmax;
    this->cropDistanceTable();
    this->uncacheCostData();
    this->uncacheRExtent();
}

double Crystal::getRmax() const
{
    double rv = (this->_rmax > 0) ?
        this->_rmax : this->_full_distance_table->maxDistance();
    return rv;
}


// r-range extended by circum diameter of the primitive cell
pair<double,double> Crystal::getRExtent() const
{
    if (!this->_rextent_cached)     this->updateRExtent();
    return this->_rextent;
}


// The equivalent unit cell vectors are mapped to fractional coordinates
// from 0 <= xi < 1.  Due to round-off errors, some vectors can end very close
// to 1.  Therefore the r-extent has to equal the largest diagonal in the cell.

void Crystal::updateRExtent() const
{
    const Lattice& lat = this->getLattice();
    double max_ucd = lat.ucMaxDiagonalLength();
    max_ucd = max_ucd*(1.0 + eps_distance) + eps_distance;
    // rmin is not defined
    double rmin = 0.0;
    double rextlo = rmin - max_ucd;
    double rexthi = this->getRmax() + max_ucd;
    // cache the results
    this->_rextent = make_pair(rextlo, rexthi);
    this->_rextent_cached = true;
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


void Crystal::Add(const Atom_t& a)
{
    Molecule::Add(a);
    // shift to unit cell and take care of zero round-off
    using NS_LIGA::eps_distance;
    Atom_t* pa = this->atoms.back();
    const Lattice& lat = this->getLattice();
    R3::Vector ucl, ucs0, ucs1;
    ucl = lat.ucvFractional(lat.fractional(pa->r));
    for (int i = 0; i != R3::Ndim; ++i)
    {
        ucs0 = ucs1 = ucl;
        ucs0[i] = 0.0;
        ucs1[i] = 1.0;
        bool nearzero = 
            lat.distance(ucl, ucs0) < eps_distance ||
            lat.distance(ucl, ucs1) < eps_distance;
        if (nearzero)   ucl[i] = 0.0;
    }
    pa->r = lat.cartesian(ucl);
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
    if (this->Badness() < NS_LIGA::eps_cost)    this->ResetBadness();
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
    if (this->Badness() < NS_LIGA::eps_cost)    this->ResetBadness();
    assert(this->_count_pairs >= 0);
}


const Crystal::TriangulationAnchor&
Crystal::getLineAnchor(const RandomWeighedGenerator& rwg)
{
    assert(countAtoms() >= 1);
    static TriangulationAnchor anch;
    anch.count = 2;
    anch.B0 = anyOffsetAtomSite(rwg);
    anch.B1 = anyOffsetAtomSite(rwg);
    return anch;
}


const Crystal::TriangulationAnchor&
Crystal::getPlaneAnchor(const RandomWeighedGenerator& rwg)
{
    assert(countAtoms() >= 1);
    static TriangulationAnchor anch;
    anch.count = 3;
    anch.B0 = anyOffsetAtomSite(rwg);
    anch.B1 = anyOffsetAtomSite(rwg);
    anch.B2 = anyOffsetAtomSite(rwg);
    return anch;
}


const Crystal::TriangulationAnchor&
Crystal::getPyramidAnchor(const RandomWeighedGenerator& rwg)
{
    return Crystal::getPlaneAnchor(rwg);
}


void Crystal::resizePairMatrices(int sz)
{
    Molecule::resizePairMatrices(sz);
    size_t sznew = this->pmx_partial_costs.rows();
    this->pmx_pair_counts.resize(sznew, sznew, 0);
}

boost::python::object Crystal::newDiffPyStructure()
{
    boost::python::object stru;
    stru = this->Molecule::newDiffPyStructure();
    const Lattice& L = this->getLattice();
    boost::python::object lattice;
    lattice = stru.attr("lattice");
    boost::python::call_method<void>(lattice.ptr(), "setLatPar",
            L.a(), L.b(), L.c(), L.alpha(), L.beta(), L.gamma());
    return stru;
}

void Crystal::setFromDiffPyStructure(boost::python::object stru)
{
    namespace python = boost::python;
    boost::python::object lattice;
    lattice = stru.attr("lattice");
    double a, b, c, alpha, beta, gamma;
    a = python::extract<double>(lattice.attr("a"));
    b = python::extract<double>(lattice.attr("b"));
    c = python::extract<double>(lattice.attr("c"));
    alpha = python::extract<double>(lattice.attr("alpha"));
    beta = python::extract<double>(lattice.attr("beta"));
    gamma = python::extract<double>(lattice.attr("gamma"));
    Lattice L(a, b, c, alpha, beta, gamma);
    this->setLattice(L);
    this->setMaxAtomCount(python::len(stru));
    this->Molecule::setFromDiffPyStructure(stru);
}


// private class methods


boost::shared_ptr<Lattice> Crystal::getDefaultLattice()
{
    static boost::shared_ptr<Lattice> default_lattice;
    if (!default_lattice.get())
    {
        default_lattice.reset(new Lattice());
    }
    return default_lattice;
}


// private methods


void Crystal::init()
{
    this->_rmin = 0.0;
    this->_rmax = 0.0;
    this->_lattice = Crystal::getDefaultLattice();
    this->_full_distance_table = this->Molecule::_distance_table;
    // the default _distance_table should exist and be blank
    assert(this->_full_distance_table.get());
    this->uncacheCostData();
    this->uncacheRExtent();
    this->setDistReuse(true);
}


void Crystal::uncacheCostData() const
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


void Crystal::uncacheRExtent() const
{
    this->_rextent_cached = false;
}


void Crystal::cropDistanceTable()
{
    DistanceTable::const_iterator lo, hi;
    lo = this->_full_distance_table->begin();
    hi = upper_bound(this->_full_distance_table->begin(),
            this->_full_distance_table->end(), this->getRmax());
    vector<double> vcropped(lo, hi);
    // Molecule::setDistanceTable must be passed DistanceTable
    // otherwise there would be infinite loop.
    DistanceTable dcropped(vcropped);
    this->Molecule::setDistanceTable(dcropped);
}


// position of random atom site offset by random lattice vector
const R3::Vector&
Crystal::anyOffsetAtomSite(const RandomWeighedGenerator& rwg) const
{
    assert(countAtoms() >= 1);
    static R3::Vector rv;
    int idx = rwg.weighedInt();
    R3::Vector mno(randomInt(2), randomInt(2), randomInt(2));
    rv = this->atoms[idx]->r + getLattice().cartesian(mno);
    return rv;
}


// End of file
