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

    // line distance table
    double pdLine[1] = {1.0};
    SandSphere ssLine(100, 1, pdLine);
    ssLine.SetGridTol(1.1);
    // unit line coordinates
    double px0[2] = {-.5,  .5};
    double py0[2] = { .0,  .0};
    double pz0[2] = { .0,  .0};
    Molecule line(&ssLine, 2, px0, py0, pz0);
    // evaluate and print fitness for line
    cout << "line:\n";
    cout << line;
    line.PrintBadness();
    line.PrintFitness();
    cout << endl;

    // square distance table
    double pdSquare[6] = {1.0, 1.0, 1.4142, 1.4142, 1.0, 1.0};
    SandSphere ssSquare(100, 6, pdSquare);
    ssSquare.SetGridTol(1.1);
    // build correct square molecule
    double px1[4] = {-.5,  .5,  .5, -.5};
    double py1[4] = {-.5, -.5,  .5,  .5};
    double pz1[4] = { .0,  .0,  .0,  .0};
    Molecule square1(&ssSquare, 4, px1, py1, pz1);
    // build defective square molecule
    double px2[4] = {-.5,  .5,  .5, -.5};
    double py2[4] = {-.5, -.5,  .5,  .0};
    double pz2[4] = { .0,  .0,  .0,  .0};
    Molecule square2(&ssSquare, 4, px2, py2, pz2);
    // evaluate and print fitness for square1
    cout << "square1:\n";
    cout << square1;
    square1.PrintBadness();
    square1.PrintFitness();
    cout << endl;
    // evaluate and print fitness for square2
    cout << "square2:\n";
    cout << square2;
    square2.PrintBadness();
    square2.PrintFitness();
    cout << endl;

    return EXIT_SUCCESS;
}
