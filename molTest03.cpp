/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test Molecule reading and writing functions
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

    Molecule mol1(&ssHex), mol2(&ssHex);
    mol1.Clear();
    mol1.ReadXY("hexagon.dxy");

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
    return EXIT_SUCCESS;
}
