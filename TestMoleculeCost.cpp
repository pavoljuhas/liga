/***********************************************************************
* Short Title: unit tests of const evaluation in Molecule
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <stdexcept>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "Molecule.hpp"

using namespace std;

class TestMoleculeCost : public CppUnit::TestFixture
{

    CPPUNIT_TEST_SUITE(TestMoleculeCost);
    CPPUNIT_TEST(test_line);
    CPPUNIT_TEST(test_square);
    CPPUNIT_TEST(test_bad_square);
    CPPUNIT_TEST_SUITE_END();

private:

    double double_eps;
    DistanceTable dst_line;
    DistanceTable dst_square;

public:

    void setUp()
    {
	double_eps = 1.0e-6;
	double line_data[1] = { 1.0 };
	dst_line = DistanceTable(line_data, 1);
	double square_data[6] = { 1.0, 1.0, 1.0, 1.0, sqrt(2.0), sqrt(2.0) };
	dst_square = DistanceTable(square_data, 6);
    }

    void test_line()
    {
	Molecule line(dst_line);
	line.Add(-0.5, 0.0, 0.0);
	line.Add(+0.5, 0.0, 0.0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, line.NormBadness(), double_eps);
    }

    void test_square()
    {
	Molecule square(dst_square);
	square.Add(-0.5, -0.5, 0.0);
	square.Add(+0.5, -0.5, 0.0);
	square.Add(+0.5, +0.5, 0.0);
	square.Add(-0.5, +0.5, 0.0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, square.NormBadness(), double_eps);
    }

    void test_bad_square()
    {
	Molecule bad_square(dst_square);
	bad_square.Add(-0.5, -0.5, 0.0);
	bad_square.Add(+0.5, -0.5, 0.0);
	bad_square.Add(+0.5, +0.5, 0.0);
	bad_square.Add(+0.0, +0.0, 0.0);
	double expectcost[4] = {
	    pow(1.0 - sqrt(0.5), 2)/2,
	    pow(1.0 - sqrt(0.5), 2)/2,
	    pow(sqrt(2.0) - sqrt(0.5), 2)/2,
	    0.0
	};
	expectcost[3] = accumulate(expectcost, expectcost + 3, 0.0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(expectcost[0],
		bad_square.getAtom(0).Badness(), double_eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(expectcost[1],
		bad_square.getAtom(1).Badness(), double_eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(expectcost[2],
		bad_square.getAtom(2).Badness(), double_eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(expectcost[3],
		bad_square.getAtom(3).Badness(), double_eps);
	double totalcost = accumulate(expectcost, expectcost + 4, 0.0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(totalcost,
		bad_square.Badness(), double_eps);
    }

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TestMoleculeCost);

// End of file
