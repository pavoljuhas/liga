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
* class TestAtomOverlapCost
*
* Comments: unit tests for AtomOverlapCost class
*
* $Id$
*
***********************************************************************/

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "AtomOverlapCost.hpp"
#include "DistanceTable.hpp"
#include "Molecule.hpp"

using namespace std;

class TestAtomOverlapCost : public CppUnit::TestFixture
{

    CPPUNIT_TEST_SUITE(TestAtomOverlapCost);
    CPPUNIT_TEST(test_notouch);
    CPPUNIT_TEST(test_overlap);
    CPPUNIT_TEST_SUITE_END();

private:

    double double_eps;
    DistanceTable dst;
    Molecule mnotouch;
    Molecule moverlap;
    AtomOverlapCost aoc;

public:

    TestAtomOverlapCost() : CppUnit::TestFixture(), aoc(&mnotouch)
    { }

    void setUp()
    {
	double dst_data[1] = { 1.0 };
	dst = DistanceTable(dst_data, 1);
        mnotouch.setDistanceTable(dst);
        mnotouch.Clear();
        mnotouch.AddAt("C", -0.5, 0.0, 0.0);
        mnotouch.AddAt("C", +0.5, 0.0, 0.0);
        moverlap = mnotouch;
        // nasty hack while there is no way to set atom radii
        Atom_t& a0 = const_cast<Atom_t&>(moverlap.getAtom(0));
        a0.radius = 2.0;
        Atom_t& a1 = const_cast<Atom_t&>(moverlap.getAtom(1));
        a1.radius = 2.0;
    }


    void test_notouch()
    {
        const Atom_t& a1 = mnotouch.getAtom(1);
        mnotouch.Pop(1);
        aoc.resetFor(&mnotouch);
        aoc.eval(&a1);
        CPPUNIT_ASSERT_EQUAL(0.0, aoc.totalCost());
        aoc.eval(&a1, AtomCost::GRADIENT);
        R3::Vector g = aoc.gradient();
        CPPUNIT_ASSERT_EQUAL(0.0, R3::norm(g));
    }


    void test_overlap()
    {
        Atom_t a1 = moverlap.getAtom(1);
        moverlap.Pop(1);
        // check if radii are properly set
        CPPUNIT_ASSERT_EQUAL(2.0, a1.radius);
        CPPUNIT_ASSERT_EQUAL(2.0, moverlap.getAtom(0).radius);
        aoc.resetFor(&moverlap);
        aoc.eval(a1);
        CPPUNIT_ASSERT_EQUAL(9.0, aoc.totalCost());
        aoc.eval(a1, AtomCost::GRADIENT);
        R3::Vector g;
        g = aoc.gradient();
        CPPUNIT_ASSERT_EQUAL(-6.0, g[0]);
        CPPUNIT_ASSERT_EQUAL(0.0, g[1]);
        CPPUNIT_ASSERT_EQUAL(0.0, g[2]);
        a1.r = -0.5, -1, 0;
        aoc.eval(a1, AtomCost::GRADIENT);
        CPPUNIT_ASSERT_EQUAL(9.0, aoc.totalCost());
        g = aoc.gradient();
        CPPUNIT_ASSERT_EQUAL(0.0, g[0]);
        CPPUNIT_ASSERT_EQUAL(6.0, g[1]);
        CPPUNIT_ASSERT_EQUAL(0.0, g[2]);
    }

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TestAtomOverlapCost);

// End of file
