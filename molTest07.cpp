/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test of Molecule::Rotate()
*
* $Id$
***********************************************************************/

#include <unistd.h>
#include "BGAlib.hpp"

int main(int argc, char *argv[])
{
    using namespace std;
    SandSphere ss(100, "square.dss");
    ss.SetGridTol(1.05);
    Molecule mol(&ss); mol.ReadXY("square.dxy");

    cout << "original molecule:" << endl;
    cout << mol; mol.PrintBadness(); cout << endl;

    cout << "applying Rotate(M_PI/4):" << endl;
    mol.Rotate(M_PI/4);
    cout << mol; mol.PrintBadness(); cout << endl;

    cout << "applying Rotate(M_PI/4, 0, 49):" << endl;
    mol.Rotate(M_PI/4, 0, 49);
    cout << mol; mol.PrintBadness(); cout << endl;

    return EXIT_SUCCESS;
}
