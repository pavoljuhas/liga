/***********************************************************************
*
* Liga Algorithm    for structure determination from pair distances
*                   Pavol Juhas
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* class TestAtomOverlapCostCrystal
*
* Comments: unit tests for AtomOverlapCostCrystal class
*
* $Id$
*
***********************************************************************/

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "AtomOverlapCostCrystal.hpp"
#include "DistanceTable.hpp"
#include "Lattice.hpp"
#include "Crystal.hpp"
#include "LigaUtils.hpp"

using namespace std;

class TestAtomOverlapCostCrystal : public CppUnit::TestFixture
{

    CPPUNIT_TEST_SUITE(TestAtomOverlapCostCrystal);
    CPPUNIT_TEST(test_notouch);
    CPPUNIT_TEST(test_overlap);
    CPPUNIT_TEST_SUITE_END();

private:

    double eps;
    DistanceTable dst;
    Lattice cubic;
    Crystal crnotouch;
    Crystal croverlap;
    AtomOverlapCostCrystal aoc;

public:

    TestAtomOverlapCostCrystal() : CppUnit::TestFixture(),
        cubic(1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
        aoc(&crnotouch)
    { }

    void setUp()
    {
        // round-off tolerance
        eps = sqrt(DOUBLE_EPS);
	double dst_data[1] = { 1.0 };
	dst = DistanceTable(dst_data, 1);
        crnotouch.setDistanceTable(dst);
        crnotouch.setLattice(cubic);
        crnotouch.setChemicalFormula("C2");
        crnotouch.Clear();
        crnotouch.AddAt("C", 0.0, 0.0, 0.0);
        crnotouch.AddAt("C", 0.5, 0.5, 0.5);
        croverlap = crnotouch;
        // nasty hack while there is no way to set atom radii
        Atom_t& a0 = const_cast<Atom_t&>(croverlap.getAtom(0));
        a0.radius = 0.5;
        Atom_t& a1 = const_cast<Atom_t&>(croverlap.getAtom(1));
        a1.radius = 0.5;
    }


    void test_notouch()
    {
        Atom_t a1 = crnotouch.getAtom(1);
        crnotouch.Pop(1);
        aoc.resetFor(&crnotouch);
        aoc.eval(a1);
        CPPUNIT_ASSERT_EQUAL(0.0, aoc.totalCost());
        aoc.eval(&a1, AtomCost::GRADIENT);
        R3::Vector g = aoc.gradient();
        CPPUNIT_ASSERT_EQUAL(0.0, R3::norm(g));
    }


    void test_overlap()
    {
        Atom_t a1 = croverlap.getAtom(1);
        croverlap.Pop(1);
        // check if radii are properly set
        CPPUNIT_ASSERT_EQUAL(0.5, a1.radius);
        CPPUNIT_ASSERT_EQUAL(0.5, croverlap.getAtom(0).radius);
        aoc.resetFor(&croverlap);
        aoc.eval(a1);
        double c;
        c = 8 * (1.75 - sqrt(3.0));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(c, aoc.totalCost(), eps);
        // gradient is zero due to bcc symmetry
        aoc.eval(a1, AtomCost::GRADIENT);
        R3::Vector g;
        g = aoc.gradient();
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, R3::norm(g), eps);
        // move it so it only touches one ball
        a1.r = -0.25, 0.0, 0.0;
        aoc.eval(a1, AtomCost::GRADIENT);
        c = 0.25 * 0.25 + 0.75 * 0.75;
        CPPUNIT_ASSERT_EQUAL(c, aoc.totalCost());
        g = aoc.gradient();
        CPPUNIT_ASSERT_EQUAL(1.0, g[0]);
        CPPUNIT_ASSERT_EQUAL(0.0, g[1]);
        CPPUNIT_ASSERT_EQUAL(0.0, g[2]);
        // check sign switching
        a1.r = +0.25, 0.0, 0.0;
        aoc.eval(a1, AtomCost::GRADIENT);
        CPPUNIT_ASSERT_EQUAL(c, aoc.totalCost());
        g = aoc.gradient();
        CPPUNIT_ASSERT_EQUAL(-1.0, g[0]);
        CPPUNIT_ASSERT_EQUAL(0.0, g[1]);
        CPPUNIT_ASSERT_EQUAL(0.0, g[2]);
    }

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TestAtomOverlapCostCrystal);

// End of file
