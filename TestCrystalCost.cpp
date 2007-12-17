/***********************************************************************
* Short Title: unit tests of cost evaluation in Crystal
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

#include "Lattice.hpp"
#include "Crystal.hpp"

using namespace std;

class TestCrystalCost : public CppUnit::TestFixture
{

    CPPUNIT_TEST_SUITE(TestCrystalCost);
    CPPUNIT_TEST(test_cube);
    CPPUNIT_TEST(test_bcc);
    CPPUNIT_TEST(test_fcc);
    CPPUNIT_TEST(test_fcc_rhomb);
    CPPUNIT_TEST_SUITE_END();

private:

    Crystal crst;
    auto_ptr<Lattice> cubic;
    auto_ptr<Lattice> rhombohedral;
    vector<double> dst_cube;
    vector<double> dst_bcc;
    vector<double> dst_fcc;

public:

    void setUp()
    {
        crst.Clear();
        crst.setRRange(0.0, 3.05);
        crst.setMaxNAtoms(4);
        // do the rest only once
        cubic.reset(new Lattice(1, 1, 1, 90, 90, 90));
        rhombohedral.reset(new Lattice(1, 1, 1, 60, 60, 60));
        // distance data up to first distance over 3
        double cube_data[] = { 0.0, 1.0, 1.41421356237, 1.73205080757, 2.0,
            2.2360679775, 2.44948974278, 2.82842712475, 3.0, 3.16227766017 };
        double bcc_data[] = { 0.0, 0.866025403784, 1.0, 1.41421356237,
            1.65831239518, 1.73205080757, 2.0, 2.17944947177, 2.2360679775,
            2.44948974278, 2.59807621135, 2.82842712475, 2.95803989155, 3.0,
            3.16227766017 };
        double fcc_data[] = { 0.0, 0.707106781187, 1.0, 1.22474487139,
            1.41421356237, 1.58113883008, 1.73205080757, 1.87082869339, 2.0,
            2.12132034356, 2.2360679775, 2.34520787991, 2.44948974278,
            2.5495097568, 2.73861278753, 2.82842712475, 2.91547594742, 3.0,
            3.08220700148 };
        size_t len_cube_data = sizeof(cube_data)/sizeof(double);
        size_t len_bcc_data = sizeof(bcc_data)/sizeof(double);
        size_t len_fcc_data = sizeof(fcc_data)/sizeof(double);
        dst_cube.assign(cube_data, cube_data + len_cube_data);
        dst_bcc.assign(bcc_data, bcc_data + len_bcc_data);
        dst_fcc.assign(fcc_data, fcc_data + len_fcc_data);
    }

    void test_cube()
    {
        crst.setLattice(*cubic);
        crst.setDistanceTable(dst_cube);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.0, 0.0, 0.0);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.1, 0.2, 0.3);
        CPPUNIT_ASSERT(crst.NormBadness() > 0.0);
        crst.Pop(1);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
    }

    void test_bcc()
    {
        crst.setLattice(*cubic);
        crst.setDistanceTable(dst_bcc);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.0, 0.0, 0.0);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.1, 0.2, 0.3);
        CPPUNIT_ASSERT(crst.NormBadness() > 0.0);
        crst.Pop(1);
        crst.Add(0.5, 0.5, 0.5);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
    }

    void test_fcc()
    {
        crst.setLattice(*cubic);
        crst.setDistanceTable(dst_fcc);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.0, 0.0, 0.0);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.5, 0.5, 0.5);
        CPPUNIT_ASSERT(crst.NormBadness() > 0.0);
        crst.Pop(1);
        crst.Add(0.5, 0.5, 0.0);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.5, 0.0, 0.5);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.0, 0.5, 0.5);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.setLattice(*rhombohedral);
        CPPUNIT_ASSERT(crst.NormBadness() > 0.0);
        crst.setLattice(*cubic);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
    }

    void test_fcc_rhomb()
    {
        Crystal crst;
        crst.setLattice(*rhombohedral);
        crst.setDistanceTable(dst_fcc);
        crst.setRRange(0.0, 3.05);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.0, 0.0, 0.0);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
        crst.Add(0.1, 0.2, 0.3);
        CPPUNIT_ASSERT(crst.NormBadness() > 0.0);
        crst.Pop(1);
        CPPUNIT_ASSERT_EQUAL(0.0, crst.NormBadness());
    }

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TestCrystalCost);

// End of file
