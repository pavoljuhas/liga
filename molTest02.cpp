/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: read distance table and coordinates for 100 random points
*     inside R=1 sphere, calculate molecular and atomic fitness/badness
*     of the ideal and disrupted molecule for a set of grid tolerances
*
* $Id$
***********************************************************************/

#include "BGAlib.hpp"

int main()
{
    using namespace std;

    double GridTolerances[] = { 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
	1.1, 1.2, 1.3, 1.4, 1.5};
    int NGridTolerances = sizeof(GridTolerances)/sizeof(double);

    // distance table
    char *distfile = "grid125.dss";
    SandSphere ssRand(500, distfile);

    // read in coordinates of random molecule
    Molecule mol1(&ssRand);
    mol1.ReadXYZ("grid125.xyz");

    // build defective random molecule
    Molecule mol2 = mol1;
    // just move atom 0 to 0,0,0
    mol2.Pop(0).Add(0,0,0);

    for (int igt = 0; igt < NGridTolerances; ++igt)
    {
	ssRand.SetGridTol(GridTolerances[igt]);
	// evaluate and print fitness for mol1
	cout << "mol1:  GridTol() = " << ssRand.GridTol() << '\t';
	cout << "MBadness() = " << mol1.Badness() << endl;
	mol1.PrintBadness();
	mol1.PrintFitness();
	cout << endl;

	// evaluate and print fitness for mol2
	cout << "mol2:  GridTol() = " << ssRand.GridTol() << '\t';
	cout << "MBadness() = " << mol2.Badness() << endl;
	mol2.PrintBadness();
	mol2.PrintFitness();
	cout << endl;
    }

    return EXIT_SUCCESS;
}
