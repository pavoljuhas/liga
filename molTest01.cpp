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
    ssSquare.SetGridTol(1.1);

    // build correct square molecule
    double px1[4] = {-.5,  .5,  .5, -.5};
    double py1[4] = {-.5, -.5,  .5,  .5};
    double pz1[4] = { .0,  .0,  .0,  .0};
    Molecule mol1(&ssSquare, 4, px1, py1, pz1);

    // build defective square molecule
    double px2[4] = {-.5,  .5,  .5, -.5};
    double py2[4] = {-.5, -.5,  .5,  .0};
    double pz2[4] = { .0,  .0,  .0,  .0};
    Molecule mol2(&ssSquare, 4, px2, py2, pz2);

    // evaluate and print fitness for mol1
    cout << "mol1:\n";
    cout << "Badness() = " << mol1.Badness() << endl;
    cout << "Fitness() = " << mol1.Fitness() << endl;
    mol1.PrintBadness();
    mol1.PrintFitness();
    cout << endl;

    // evaluate and print fitness for mol2
    cout << "mol2:\n";
    cout << "Badness() = " << mol2.Badness() << endl;
    cout << "Fitness() = " << mol2.Fitness() << endl;
    mol2.PrintBadness();
    mol2.PrintFitness();
    cout << endl;

    return EXIT_SUCCESS;
}
