/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test Molecule reading and writing functions
*
* $Id$
***********************************************************************/

#include "BGAlib.hpp"

int main()
{
    using namespace std;

    // distance table
    SandSphere ssHex(100, "hexagon.dss");

    Molecule mol(&ssHex);
    mol.Clear();
    mol.ReadXY("hexagon.dxy");

    cout << "mol:  MBadness() = " << mol.MBadness() << endl;
    for (int i = 0; i < mol.NAtoms; ++i)
    {
	cout << "ABadness(" << i << ") = " << mol.ABadness(i) << "  ";
	cout << "AFitness(" << i << ") = " << mol.AFitness(i) << '\n';
    }
    cout << endl;

    cout << "mol.OutFmtGrid():" << endl;
    cout << mol.OutFmtGrid() << endl;

    cout << "mol.OutFmtXY():" << endl;
    cout << mol.OutFmtXY() << endl;

    cout << "mol.OutFmtAeye():" << endl;
    cout << mol.OutFmtAeye() << endl;

    return EXIT_SUCCESS;
}
