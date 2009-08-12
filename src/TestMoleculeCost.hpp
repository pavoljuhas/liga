/***********************************************************************
* Short Title: unit tests of const evaluation in Molecule
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <cxxtest/TestSuite.h>

#include "Molecule.hpp"

using namespace std;

class TestMoleculeCost : public CxxTest::TestSuite
{
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
            Molecule line;
            line.setDistanceTable(dst_line);
            line.AddAt("", -0.5, 0.0, 0.0);
            line.AddAt("", +0.5, 0.0, 0.0);
            TS_ASSERT_DELTA(0.0, line.cost(), double_eps);
        }


        void test_square()
        {
            Molecule square;
            square.setDistanceTable(dst_square);
            square.AddAt("", -0.5, -0.5, 0.0);
            square.AddAt("", +0.5, -0.5, 0.0);
            square.AddAt("", +0.5, +0.5, 0.0);
            square.AddAt("", -0.5, +0.5, 0.0);
            TS_ASSERT_DELTA(0.0, square.cost(), double_eps);
        }


        void test_bad_square()
        {
            Molecule bad_square;
            bad_square.setDistanceTable(dst_square);
            bad_square.AddAt("", -0.5, -0.5, 0.0);
            bad_square.AddAt("", +0.5, -0.5, 0.0);
            bad_square.AddAt("", +0.5, +0.5, 0.0);
            bad_square.AddAt("", +0.0, +0.0, 0.0);
            double expectcost[4] = {
                pow(1.0 - sqrt(0.5), 2)/2,
                pow(1.0 - sqrt(0.5), 2)/2,
                pow(sqrt(2.0) - sqrt(0.5), 2)/2,
                0.0
            };
            expectcost[3] = accumulate(expectcost, expectcost + 3, 0.0);
            TS_ASSERT_DELTA(expectcost[0],
                    bad_square.getAtom(0).Badness(), double_eps);
            TS_ASSERT_DELTA(expectcost[1],
                    bad_square.getAtom(1).Badness(), double_eps);
            TS_ASSERT_DELTA(expectcost[2],
                    bad_square.getAtom(2).Badness(), double_eps);
            TS_ASSERT_DELTA(expectcost[3],
                    bad_square.getAtom(3).Badness(), double_eps);
            double totalcost = accumulate(expectcost, expectcost + 4, 0.0);
            TS_ASSERT_DELTA(totalcost, bad_square.Badness(), double_eps);
            // make sure recalculate returns the same non-zero cost
            bad_square.recalculate();
            TS_ASSERT_DELTA(totalcost, bad_square.Badness(), double_eps);
        }

};  // class TestMoleculeCost

// End of file
