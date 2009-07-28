/***********************************************************************
* Short Title: unit tests of least squares atom relaxation
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

#include "AtomCost.hpp"
#include "LigaUtils.hpp"
#include "Crystal.hpp"

using namespace std;

class TestRelaxCrystalAtom : public CppUnit::TestFixture
{

    CPPUNIT_TEST_SUITE(TestRelaxCrystalAtom);
    CPPUNIT_TEST(test_RelaxBCC);
    CPPUNIT_TEST(test_RelaxFCC);
    CPPUNIT_TEST(test_RelaxTriclinic);
    CPPUNIT_TEST_SUITE_END();

private:

    double eps_distance;
    auto_ptr<DistanceTable> dtgt_bcc;
    auto_ptr<DistanceTable> dtgt_fcc;
    auto_ptr<DistanceTable> dtgt_triclinic;
    auto_ptr<Crystal> crbcc;
    auto_ptr<Crystal> crfcc;
    auto_ptr<Crystal> crtriclinic;

    template <class T> T* newFromFile(const string& filename)
    {
        T instance;
        ifstream fid(filename.c_str());
        fid >> instance;
        T* rv = new T(instance); 
        return rv;
    }

public:

    void setUp()
    {
	eps_distance = NS_LIGA::eps_distance;
        // load distance tables
        DistanceTable* ptbl;
        ptbl = this->newFromFile<DistanceTable>("crystals/bcc.dst");
        this->dtgt_bcc.reset(ptbl);
        ptbl = this->newFromFile<DistanceTable>("crystals/fcc.dst");
        this->dtgt_fcc.reset(ptbl);
        ptbl = this->newFromFile<DistanceTable>("crystals/triclinic.dst");
        this->dtgt_triclinic.reset(ptbl);
        // bcc
        this->crbcc.reset(new Crystal(*this->dtgt_bcc));
        this->crbcc->ReadFile("crystals/bcc.stru");
        // fcc
        this->crfcc.reset(new Crystal(*this->dtgt_fcc));
        this->crfcc->ReadFile("crystals/fcc.stru");
        // triclinic
        this->crtriclinic.reset(new Crystal(*this->dtgt_triclinic));
        this->crtriclinic->ReadFile("crystals/triclinic.stru");
    }

    void tearDown() 
    { }

    void test_RelaxBCC()
    {
        Atom_t a1 = this->crbcc->getAtom(1);
        this->crbcc->Pop(1);
        Atom_t arx = a1;
        R3::Vector offset(0.013, -0.07, -0.03);
        arx.r += offset;
        this->crbcc->RelaxExternalAtom(arx);
        double drx = R3::distance(a1.r, arx.r);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, drx, eps_distance);
    }

    void test_RelaxFCC()
    {
        Atom_t a1 = this->crfcc->getAtom(1);
        this->crfcc->Pop(1);
        Atom_t arx = a1;
        R3::Vector offset(0.013, -0.07, -0.03);
        arx.r += offset;
        this->crfcc->RelaxExternalAtom(arx);
        double drx = R3::distance(a1.r, arx.r);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, drx, eps_distance);
    }

    void test_RelaxTriclinic()
    {
        Atom_t a1 = this->crtriclinic->getAtom(1);
        this->crtriclinic->Pop(1);
        Atom_t arx = a1;
        R3::Vector offset(0.0013, -0.007, -0.003);
        arx.r += offset;
        this->crtriclinic->RelaxExternalAtom(arx);
        double drx = R3::distance(a1.r, arx.r);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, drx, eps_distance);
    }


};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TestRelaxCrystalAtom);

// End of file
