/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: create correct and defect square molecule and evaluate
*     molecular and atomic fitness/badness
*
* $Id$
***********************************************************************/

#include "BGAlib.hpp"

int main()
{
    using namespace std;
    // square distance table
    double pdSquare[6] = {1.0, 1.0, 1.4142, 1.4142, 1.0, 1.0};
    SandSphere ssSquare(100, 6, pdSquare);

    // build correct square molecule
    double px1[4] = {-.5,  .5,  .5, -.5};
    double py1[4] = {-.5, -.5,  .5,  .5};
    Molecule mol1(&ssSquare, 4, px1, py1);

    // build defective square molecule
    double px2[4] = {-.5,  .5,  .5, -.5};
    double py2[4] = {-.5, -.5,  .5,  .5};
    Molecule mol2(&ssSquare, 4, px2, py2);

    // evaluate and print fitness for mol1
    cout << "mol1:\n";
    cout << "MBadness() = " << mol1.MBadness() << endl;
    cout << "MFitness() = " << mol1.MFitness() << endl;
    for (int i = 0; i < mol1.NAtoms; ++i)
    {
	cout << "ABadness(" << i << ") = " << mol1.ABadness(i) << '\t';
	cout << "AFitness(" << i << ") = " << mol1.AFitness(i) << '\n';
    }
    cout << endl;

    // evaluate and print fitness for mol2
    cout << "mol2:\n";
    cout << "MBadness() = " << mol2.MBadness() << endl;
    cout << "MFitness() = " << mol2.MFitness() << endl;
    for (int i = 0; i < mol2.NAtoms; ++i)
    {
	cout << "ABadness(" << i << ") = " << mol2.ABadness(i) << '\t';
	cout << "AFitness(" << i << ") = " << mol2.AFitness(i) << '\n';
    }
    cout << endl;

    return EXIT_SUCCESS;
}
