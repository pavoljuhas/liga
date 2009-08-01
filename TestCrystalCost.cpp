/***********************************************************************
* Short Title: unit tests of cost evaluation in Crystal
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "Lattice.hpp"
#include "Crystal.hpp"
#include "LigaUtils.hpp"
#include "AtomCost.hpp"

using namespace std;

class TestCrystalCost : public CppUnit::TestFixture
{

    CPPUNIT_TEST_SUITE(TestCrystalCost);
    CPPUNIT_TEST(test_cube);
    CPPUNIT_TEST(test_bcc);
    CPPUNIT_TEST(test_fcc);
    CPPUNIT_TEST(test_fcc_rhomb);
    CPPUNIT_TEST(test_gradient_bcc);
    CPPUNIT_TEST(test_gradient_triclinic);
    CPPUNIT_TEST_SUITE_END();

private:

    double double_eps;
    double gradient_eps;
    Crystal crst;
    auto_ptr<Lattice> cubic;
    auto_ptr<Lattice> rhombohedral;
    vector<double> dst_cube;
    vector<double> dst_bcc;
    vector<double> dst_fcc;

public:

    void setUp()
    {
	double_eps = DOUBLE_EPS;
	gradient_eps = 1e-6;
        crst.Clear();
        crst.setRmax(3.05);
        ChemicalFormula::value_type elcnt("C", 4);
        ChemicalFormula formula(1, elcnt);
        crst.setChemicalFormula(formula);
        // do the rest only once
        if (!cubic.get())   cubic.reset(new Lattice(1, 1, 1, 90, 90, 90));
        if (!rhombohedral.get())    rhombohedral.reset(
                new Lattice(sqrt(0.5), sqrt(0.5), sqrt(0.5), 60, 60, 60));
        // distance data up to first distance over 3
        double cube_data[] = { 1.0, 1.41421356237, 1.73205080757, 2.0,
            2.2360679775, 2.44948974278, 2.82842712475, 3.0, 3.16227766017 };
        double bcc_data[] = { 0.866025403784, 1.0, 1.41421356237,
            1.65831239518, 1.73205080757, 2.0, 2.17944947177, 2.2360679775,
            2.44948974278, 2.59807621135, 2.82842712475, 2.95803989155, 3.0,
            3.16227766017 };
        double fcc_data[] = { 0.707106781187, 1.0, 1.22474487139,
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
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.0, 0.0, 0.0);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.1, 0.2, 0.3);
        CPPUNIT_ASSERT(crst.cost() > 0.0);
        crst.Pop(1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
    }

    void test_bcc()
    {
        crst.setLattice(*cubic);
        crst.setDistanceTable(dst_bcc);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.0, 0.0, 0.0);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.1, 0.2, 0.3);
        CPPUNIT_ASSERT(crst.cost() > 0.0);
        crst.Pop(1);
        crst.AddAt("C", 0.5, 0.5, 0.5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
    }

    void test_fcc()
    {
        crst.setLattice(*cubic);
        crst.setDistanceTable(dst_fcc);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.0, 0.0, 0.0);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.5, 0.5, 0.5);
        CPPUNIT_ASSERT(crst.cost() > 0.0);
        crst.Pop(1);
        crst.AddAt("C", 0.5, 0.5, 0.0);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.5, 0.0, 0.5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.0, 0.5, 0.5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.setLattice(*rhombohedral);
        CPPUNIT_ASSERT(crst.cost() > 0.0);
        crst.setLattice(*cubic);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
    }

    void test_fcc_rhomb()
    {
        crst.setLattice(*rhombohedral);
        crst.setDistanceTable(dst_fcc);
        crst.setRmax(3.05);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.0, 0.0, 0.0);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
        crst.AddAt("C", 0.1, 0.2, 0.3);
        CPPUNIT_ASSERT(crst.cost() > 0.0);
        crst.Pop(1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, crst.cost(), double_eps);
    }

    R3::Vector analytical_gradient(const Atom_t& a0)
    {
        AtomCost* atomcost = crst.getAtomCostCalculator();
        atomcost->eval(a0, AtomCost::GRADIENT);
        R3::Vector grad = atomcost->gradient();
        return grad;
    }

    R3::Vector numerical_gradient(const Atom_t& a0)
    {
        const double delta = 1e-8;
        AtomCost* atomcost = crst.getAtomCostCalculator();
        double ac0 = atomcost->eval(a0);
        R3::Vector numgrad;
        for (int i = 0; i < 3; ++i)
        {
            Atom_t a1 = a0;
            a1.r[i] += delta;
            double ac1 = atomcost->eval(a1);
            numgrad[i] = (ac1 - ac0) / delta;
        }
        return numgrad;
    }

    void test_gradient_bcc()
    {
        crst.setLattice(*cubic);
        crst.setDistanceTable(dst_bcc);
        crst.AddAt("C", 0.0, 0.0, 0.0);
        Atom_t a1("C", 0.5, 0.5, 0.5);
        R3::Vector ga = analytical_gradient(a1);
        R3::Vector gn = numerical_gradient(a1);
        double gdiff = R3::distance(gn, ga);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, gdiff, gradient_eps);
    }

    void test_gradient_triclinic()
    {
        DistanceTable dst;
        ifstream fid("crystals/triclinic.dst");
        fid >> dst;
        crst.setDistanceTable(dst);
        crst.ReadFile("crystals/triclinic.stru");
        Atom_t a1 = crst.getAtom(1);
        crst.Pop(1);
        R3::Vector ga = analytical_gradient(a1);
        double ganorm = R3::norm(ga);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, ganorm, gradient_eps);
        Atom_t a2("C", 0.5, 0.5, 0.5);
        R3::Vector gn;
        ga = analytical_gradient(a2);
        gn = numerical_gradient(a2);
        double gdiff = R3::distance(gn, ga);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, gdiff, gradient_eps);
    }

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TestCrystalCost);

// End of file
