/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test Molecule operators: operator=, copy constructor,
*     Shift(), Center(), Part(), Pop(), Add(), MoveAtomTo()
*
* $Id$
***********************************************************************/

#include <unistd.h>
#include "BGAlib.hpp"

int main()
{
    using namespace std;

    // distance table
    SandSphere ssHex(100, "hexagon.dss");

    Molecule mol1(&ssHex);
    mol1.ReadXY("hexagon.dxy");

    cout << "mol1:\n";
    mol1.PrintBadness();
    cout << endl;

    Molecule mol2(&ssHex);

    cout << "mol2: testing mol2 = mol1:" << endl;
    mol2 = mol1;
    cout << mol2; mol2.PrintBadness(); cout << endl;

    cout << "testing mol2.Shift(502, 0):" << endl;
    mol2.Shift(502, 0);
    cout << mol2; mol2.PrintBadness(); cout << endl;

    cout << "testing mol2.Center():" << endl;
    mol2.Center();
    cout << mol2; mol2.PrintBadness(); cout << endl;

    cout << "testing mol2.Pop(0):" << endl;
    mol2.Pop(0);
    cout << mol2; mol2.PrintBadness(); cout << endl;

    cout << "testing mol2.Add(50, 0):" << endl;
    mol2.Add(50, 0);
    cout << mol2; mol2.PrintBadness(); cout << endl;

    cout << "testing Molecule mol3(mol2); mol3.Part(idx); mol2.Pop(idx):\n";
    int i1[3] = {0, 1, 2}; list<int> idx(i1, i1+3);
    Molecule mol3(mol2);
    mol3.Part(idx);
    mol2.Pop(idx);
    cout << "mol2:" << endl << mol2; mol2.PrintBadness();
    cout << "mol3:" << endl << mol3; mol3.PrintBadness();
    cout << endl;

    cout << "testing mol3.Add(mol2):" << endl;
    mol3.Add(mol2);
    cout << "mol3:" << endl << mol3; mol3.PrintBadness();
    cout << endl;

    cout << "testing mol3.MoveAtomTo(0, 0, 7):" << endl;
    mol3.MoveAtomTo(0, 0, 7);
    cout << "mol3:" << endl << mol3; mol3.PrintBadness();
    cout << endl;
    cout << "testing mol3.MoveAtomTo(0, 25, 43):" << endl;
    mol3.MoveAtomTo(0, 25, 43);
    cout << "mol3:" << endl << mol3; mol3.PrintBadness();
    cout << endl;

    return EXIT_SUCCESS;
}
