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
    const double pemax = 0.95;
    const double pallway = 0.01;
    const double avgmb = 0.05;
    ////////////////////////////////////////////////////////////////////////
    if (argc == 1)
    {
	cerr << "usage: " <<
	    "molTest08 distance_file.dss [snapshot_file.dxy]" << endl;
	return EXIT_SUCCESS;
    }
    // here argc > 1
    char *distance_file = argv[1];
    SandSphere* ss;
    try
    {
	ss = new SandSphere(500, distance_file);
    }
    catch (IOError)
    {
	return EXIT_FAILURE;
    }
    ss->SetGridTol(1.1);
    char *snapshot_file = NULL;
    if (argc > 2 && strlen(argv[2]) > 0)
    {
	snapshot_file = argv[2];
	cout << "molecule snapshots go to " << snapshot_file << endl;
    }

    // set lastMBadness to a maximum double
    numeric_limits<double> double_info;
    valarray<double> lastMBadness(double_info.max(), ss->NAtoms);
    double best_largest = double_info.max();
    valarray<int> improved(1, logsize);

    Molecule mol(ss);


    int maxatoms = 0;
    bool go_all_way = false;
    for (int trial = 0; ; ++trial)
    {
	// calculate pe
	double pe, impr_rate;
	if (mol.NAtoms == ss->NAtoms || mol.MBadness() > 10*mol.NAtoms)
	{
	    pe = 0.0;
	    go_all_way = false;
	}
	else if (mol.NAtoms == 0 || go_all_way)
	{
	    pe = 1.0;
	}
	else
	{
	    impr_rate = 1.0*improved.sum()/improved.size();
	    pe = impr_rate*(pemax-pemin)+pemin;
	    if (impr_rate <= 0.2 && pallway > gsl_rng_uniform(BGA::rng))
		go_all_way = true;
	}
	cout << trial;
	if (pe > gsl_rng_uniform(BGA::rng))
	{
	    mol.Evolve();
	    cout << "  Evolve()" << "  NAtoms = " << mol.NAtoms << endl;
	}
	else
	{
	    int Npop = 0;
	    for (int i = 0; i < mol.NAtoms; ++i) Npop += (mol.ABadness(i)>0);
	    Npop = (int) ceil( (mol.MBadness()/mol.NAtoms + Npop)/2.0 );
	    Npop = min(Npop, 5);
	    if (Npop > 1)
		Npop = 1 + gsl_rng_uniform_int(BGA::rng, Npop-1);
	    mol.Degenerate(Npop);
	    cout << "  Degenerate(" << Npop << ")  NAtoms = " << mol.NAtoms <<  endl;
	}
	mol.PrintBadness();
	if (mol.NAtoms == ss->NAtoms)
	    cout << "mol.MBadness()/NAtoms = " << mol.MBadness()/ss->NAtoms << endl;
	// update lastMBadness and improved
	int ilog = trial % logsize;
	if (mol.MBadness()<=lastMBadness[mol.NAtoms-1] || mol.MBadness() == 0.0)
	{
	    if (mol.NAtoms > maxatoms)
	    {
		best_largest = lastMBadness[mol.NAtoms-1];
		maxatoms = mol.NAtoms;
	    }
	    lastMBadness[mol.NAtoms-1] = mol.MBadness();
	    improved[ilog] = 1;
	    if (snapshot_file != NULL && mol.NAtoms == maxatoms &&
		    best_largest >= mol.MBadness())
	    {
		maxatoms = mol.NAtoms;
		best_largest = mol.MBadness();
		cout << "saving best molecule" << endl;
		mol.WriteXY(snapshot_file);
	    }
	    if (mol.NAtoms == ss->NAtoms)
	    {
		if (mol.MBadness() < avgmb*ss->NAtoms)
		{
		    cout << "that is solution!" << endl;
		    break;
		}
	    }
	}
	else
	{
	    improved[ilog] = 0;
	    lastMBadness[mol.NAtoms-1] = mol.MBadness();
	}
	cout << endl;
    }
    return EXIT_SUCCESS;
}
