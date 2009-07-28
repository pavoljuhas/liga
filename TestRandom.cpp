/***********************************************************************
* Short Title: unit tests for Lattice class
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <algorithm>
#include <stdexcept>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "Random.hpp"

using namespace std;
using namespace NS_LIGA;
using R3::MatricesAlmostEqual;
using R3::VectorsAlmostEqual;

class TestRandom : public CppUnit::TestFixture
{

    CPPUNIT_TEST_SUITE(TestRandom);
    CPPUNIT_TEST(test_randomPickFew);
    CPPUNIT_TEST(test_RGW_setWeights);
    CPPUNIT_TEST(test_RGW_weighedPick);
    CPPUNIT_TEST(test_RGW_weighedInt);
    CPPUNIT_TEST_SUITE_END();

private:

    RandomWeighedGenerator rwg;
    static const double precision = 1.0e-12;

public:

    void setUp()
    {
    }

    void tearDown()
    {
    }

    void test_randomPickFew()
    {
        PickType pick, pickall;
        CPPUNIT_ASSERT_THROW(randomPickFew(50, 49), out_of_range);
        size_t cnt0 = 0;
        size_t attempts = 10000;
        for (size_t i = 0; i < attempts; ++i)
        {
            pick = randomPickFew(5, 10);
            if (pick[0] == 1)   cnt0++;
        }
        double p = 0.1;
        double avgcnt0 = attempts*p;
        double sigcnt0 = sqrt(attempts*p*(1 - p));
        CPPUNIT_ASSERT(cnt0 < avgcnt0 + 6*sigcnt0);
        CPPUNIT_ASSERT(cnt0 > avgcnt0 - 6*sigcnt0);
        // test pick uniqueness
        for (size_t i = 1; i < 1000; ++i)
        {
            pickall = randomPickFew(10, 10);
            sort(pickall.begin(), pickall.end());
            CPPUNIT_ASSERT(pickall.end() ==
                    unique(pickall.begin(), pickall.end()) );
        }
        CPPUNIT_ASSERT(randomPickFew(0, 100).empty());
        CPPUNIT_ASSERT(0 == randomPickFew(1, 1)[0]);
    }

    void test_RGW_setWeights()
    {
        double badwts[3] = {1.0, -5, 2.0};
        CPPUNIT_ASSERT_THROW(rwg.setWeights(badwts, badwts + 3), out_of_range);
    }

    void test_RGW_weighedPick()
    {
        double wts[3] = {1, 2, 4};
        rwg.setWeights(wts, wts + 3);
        size_t cntfirst0 = 0;
        size_t cntsecond1 = 0;
        size_t cntthird2 = 0;
        size_t attempts = 10000;
        for (size_t i = 0; i < attempts; ++i)
        {
            const PickType& sel = rwg.weighedPick(3);
            if (sel[0] == 0)    cntfirst0++;
            if (sel[1] == 1)    cntsecond1++;
            if (sel[2] == 2)    cntthird2++;
        }
        double pfirst0 = 1.0/7;
        double psecond1 = 1.0/7 * 1.0/3 + 4.0/7 * 2.0/3;
        double pthird2 = 1.0/7 * 1.0/3 + 2.0/7 * 1.0/5;
        double avgfirst0 = pfirst0*attempts;
        double sigfirst0 = sqrt(pfirst0*attempts*(1 - pfirst0));
        double avgsecond1 = psecond1*attempts;
        double sigsecond1 = sqrt(psecond1*attempts*(1 - psecond1));
        double avgthird2 = pthird2*attempts;
        double sigthird2 = sqrt(pthird2*attempts*(1 - pthird2));
        CPPUNIT_ASSERT(cntfirst0 < avgfirst0 + 6*sigfirst0);
        CPPUNIT_ASSERT(cntfirst0 > avgfirst0 - 6*sigfirst0);
        CPPUNIT_ASSERT(cntsecond1 < avgsecond1 + 6*sigsecond1);
        CPPUNIT_ASSERT(cntsecond1 > avgsecond1 - 6*sigsecond1);
        CPPUNIT_ASSERT(cntthird2 < avgthird2 + 6*sigthird2);
        CPPUNIT_ASSERT(cntthird2 > avgthird2 - 6*sigthird2);
    }

    void test_RGW_weighedInt()
    {
        double wts[2] = {3, 1};
        rwg.setWeights(wts, wts + 2);
        size_t cnt0 = 0;
        size_t attempts = 10000;
        for (size_t i = 0; i < attempts; ++i)
        {
            if (rwg.weighedInt() == 0)  cnt0++;
        }
        double p0 = 3.0/4;
        double avgcnt0 = p0 * attempts;
        double sigcnt0 = sqrt(p0*attempts*(1 - p0));
        CPPUNIT_ASSERT(cnt0 < avgcnt0 + 6*sigcnt0);
        CPPUNIT_ASSERT(cnt0 > avgcnt0 - 6*sigcnt0);
    }

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TestRandom);

// End of file
