/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: read distance table and coordinates for 100 random points
*     inside R=1 circle, calculate molecular and atomic fitness/badness
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
    char *distfile = "rand100.dss";
    SandSphere ssRand(100, distfile);

    // read in coordinates:
    vector<double> vx1, vy1, vx2, vy2;
    ifstream fid("rand100.dxy");
    // first line contains comment
    string dummy;
    getline(fid, dummy);
    double xi, yi;
    while ( !fid.eof() && (fid >> xi >> yi) )
    {
	vx1.push_back(xi);
	vy1.push_back(yi);
    }
    fid.close();

    // build correct random molecule
    Molecule mol1(&ssRand, vx1, vy1);

    Molecule *mpj = new Molecule(&ssRand, vx1, vy1);
    cout << "*mpj: MBadness() = " << mpj->MBadness() << endl;
    

    // build defective random molecule
    vx2 = vx1;
    vy2 = vy1;
    vx2[0] = 0.0;
    vy2[0] = 0.0;
    Molecule mol2(&ssRand, vx2, vy2);

    // test destructor
    delete mpj;

    for (int igt = 0; igt < NGridTolerances; ++igt)
    {
	ssRand.SetGridTol(GridTolerances[igt]);
    // evaluate and print fitness for mol1
    cout << "mol1:  GridTol() = " << ssRand.GridTol() << '\t';
    cout << "MBadness() = " << mol1.MBadness() << endl;
    for (int i = 0; i < mol1.NAtoms; ++i)
    {
	cout << "ABadness(" << i << ") = " << mol1.ABadness(i) << "  ";
	cout << "AFitness(" << i << ") = " << mol1.AFitness(i) << '\n';
    }
    cout << endl;

    // evaluate and print fitness for mol2
    cout << "mol2:  GridTol() = " << ssRand.GridTol() << '\t';
    cout << "MBadness() = " << mol2.MBadness() << endl;
    for (int i = 0; i < mol2.NAtoms; ++i)
    {
	cout << "ABadness(" << i << ") = " << mol2.ABadness(i) << "  ";
	cout << "AFitness(" << i << ") = " << mol2.AFitness(i) << '\n';
    }
    cout << endl;
    }

    return EXIT_SUCCESS;
}
