/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test Molecule operators Shift(), Center(), Part(), Pop(),
*     Add(), MoveAtom()
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

    Molecule mol2 = mol1;
//    Molecule mol2(&ssHex);

    mol2 = mol1;
    cout << "mol2:\n";
    mol2.PrintBadness();
    cout << endl;

    cout << mol2 << endl;
    mol2.Shift(502, 0);
    cout << mol2 << endl;
    cout << "mol2:\n";
    mol2.PrintBadness();
    cout << endl;

    cout << "mol1:\n";
    mol1.PrintBadness();
    cout << endl;

    /*
    Molecule& Shift(int dh, int dk);	// shift all atoms
    Molecule& Center();			// center w/r to the center of mass
    Molecule& Part(const list<int>& cidx);	// keep only specified atoms
    Molecule& Pop(const list<int>& cidx);	// remove specified atoms
    Molecule& Pop(const int cidx);	// remove 1 atom
    Molecule& Clear();			// remove all atoms
    Molecule& Add(Molecule& m);		// add specified molecule
    Molecule& Add(int nh, int nk);	// add single atom
    Molecule& MoveAtom(int idx, int nh, int nk);	// move 1 atom

    cout << "mol1:  MBadness() = " << mol1.MBadness() << endl;
    for (int i = 0; i < mol1.NAtoms; ++i)
    {
	cout << "ABadness(" << i << ") = " << mol1.ABadness(i) << "  ";
	cout << "AFitness(" << i << ") = " << mol1.AFitness(i) << '\n';
    }
    cout << endl;

    cout << "mol1.OutFmtGrid():" << endl;
    cout << mol1.OutFmtGrid() << endl;

    cout << "mol1.OutFmtXY():" << endl;
    cout << mol1.OutFmtXY() << endl;

    cout << "mol1.OutFmtAeye():" << endl;
    cout << mol1.OutFmtAtomEye() << endl;

    cout << "reading mol2 from Grid format" << endl;
    mol1.WriteGrid("molTest03.tmp");
    ifstream fid("molTest03.tmp");
    fid >> mol2;
    fid.close();
    cout << "mol2:  MBadness() = " << mol2.MBadness() << endl;
    for (int i = 0; i < mol2.NAtoms; ++i)
    {
	cout << "ABadness(" << i << ") = " << mol2.ABadness(i) << "  ";
	cout << "AFitness(" << i << ") = " << mol2.AFitness(i) << '\n';
    }
    cout << endl;

    cout << "reading mol2 from XY format" << endl;
    mol1.WriteXY("molTest03.tmp");
    fid.open("molTest03.tmp");
    fid >> mol2;
    fid.close();
    cout << "mol2:  MBadness() = " << mol2.MBadness() << endl;
    for (int i = 0; i < mol2.NAtoms; ++i)
    {
	cout << "ABadness(" << i << ") = " << mol2.ABadness(i) << "  ";
	cout << "AFitness(" << i << ") = " << mol2.AFitness(i) << '\n';
    }
    cout << endl;

    unlink("molTest03.tmp");
    */
    return EXIT_SUCCESS;
}
