/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test Molecule::Evolve(), Molecule::Degenerate(),
*     automatic search for good configuration
*
* $Id$
***********************************************************************/

#include <limits>
#include <unistd.h>
#include "BGAlib.hpp"

int main(int argc, char *argv[])
{
    using namespace std;

    // parameters
    const int logsize = 10;
    const double pemin = 0.25;
    const double pemax = 0.9;
    const double pallway = 0.01;
    const double avgmb = 0.05;
    ////////////////////////////////////////////////////////////////////////
    if (argc == 1)
    {
	cerr << "usage: " <<
	    "molTest08 distance_file.dss [seed] [snapshot_file.xyz]" << endl;
	return EXIT_SUCCESS;
    }
    // here argc > 1
    char *distance_file = argv[1];
    DistanceTable *dtab;
    try
    {
	dtab = new DistanceTable(distance_file);
    }
    catch (IOError)
    {
	return EXIT_FAILURE;
    }
    char *snapshot_file = NULL;
    if (argc > 2 && strlen(argv[2]) > 0)
    {
	unsigned long int seed;
	seed = atoi(argv[2]);
	cout << "setting seed to " << seed << endl;
	gsl_rng_set(BGA::rng, seed);
    }
    if (argc > 3 && strlen(argv[3]) > 0)
    {
	snapshot_file = argv[3];
	cout << "molecule snapshots go to " << snapshot_file << endl;
    }

    // set lastMBadness to a maximum double
    numeric_limits<double> double_info;
    valarray<double> lastMBadness(double_info.max(), dtab->NAtoms);
    double best_largest = double_info.max();
    valarray<int> improved(1, logsize);

    Molecule mol(*dtab);
    int fileno = 0;

    int maxatoms = 0;
    bool go_all_way = false;
    for (int trial = 0; ; ++trial)
    {
	// calculate pe
	double pe, impr_rate;
	if (mol.NAtoms() == mol.max_NAtoms() || mol.Badness() > 10*mol.NAtoms())
	{
	    pe = 0.0;
	    go_all_way = false;
	}
	else if (mol.NAtoms() == 0 || go_all_way)
	{
	    pe = 1.0;
	}
	else
	{
	    impr_rate = 1.0*improved.sum()/improved.size();
	    pe = impr_rate*(pemax-pemin)+pemin;
	    if (impr_rate >= 0.66 && pallway > gsl_rng_uniform(BGA::rng))
		go_all_way = true;
	    if (mol.Badness() == 0)
		pe = pemax;
	}
	cout << trial;
	if (pe > gsl_rng_uniform(BGA::rng))
	{
	    mol.Evolve(10,20,1000);
	    cout << "  Evolve()" << "  NAtoms = " << mol.NAtoms() << endl;
	}
	else
	{
	    int Npop = 0;
	    Npop = (int) ceil(1.0*mol.Badness() / mol.NAtoms());
	    Npop = min(Npop, 5);
	    if (Npop > 1)
		Npop = 1 + gsl_rng_uniform_int(BGA::rng, Npop-1);
	    mol.Degenerate(Npop);
	    cout << "  Degenerate(" << Npop << ")  NAtoms = " << mol.NAtoms() <<  endl;
	}
	mol.PrintBadness();
	if (mol.NAtoms() == mol.max_NAtoms())
	    cout << "mol.Badness()/NAtoms = " << mol.Badness()/mol.max_NAtoms() << endl;
	// update lastMBadness and improved
	int ilog = trial % logsize;
	if (mol.Badness()<=lastMBadness[mol.NAtoms()-1] || mol.Badness() == 0.0)
	{
	    if (mol.NAtoms() > maxatoms)
	    {
		best_largest = lastMBadness[mol.NAtoms()-1];
		maxatoms = mol.NAtoms();
	    }
	    lastMBadness[mol.NAtoms()-1] = mol.Badness();
	    improved[ilog] = 1;
	    /*
	    if (snapshot_file != NULL && mol.NAtoms() == maxatoms &&
		    best_largest >= mol.Badness())
	    {
		best_largest = mol.Badness();
		cout << "saving best molecule" << endl;
		mol.WriteXYZ(snapshot_file);
	    }
	    */
	    if (snapshot_file != NULL && improved[ilog])
	    {
		cout << "saving best molecule" << endl;
		char fname[255];
		sprintf(fname, "%s%04i", snapshot_file, fileno);
		sprintf(fname, "%s", snapshot_file);
		mol.WriteAtomEye(fname);
	    }
	    if (mol.NAtoms() == mol.max_NAtoms())
	    {
		if (mol.Badness() < avgmb*mol.max_NAtoms())
		{
		    cout << "that is solution!" << endl;
		    break;
		}
	    }
	}
	else
	{
	    improved[ilog] = 0;
	    lastMBadness[mol.NAtoms()-1] = mol.Badness();
	}
	cout << endl;
    }
    return EXIT_SUCCESS;
}
