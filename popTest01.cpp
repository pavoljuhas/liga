/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: create populations with 10 empty molecules and 10 hexagons
*
* $Id$
***********************************************************************/

#include "BGAlib.hpp"

int main()
{
    using namespace std;
    // read hexagon distance table
    SandSphere ssHex(100, "hexagon.dss");

    Molecule mol(&ssHex); mol.ReadXY("hexagon.dxy");

    Population Empties(10, &ssHex);
    Population Hexagons(10, mol);

    cout << "===== Empties =====" << endl;
    for (Population::iterator i = Empties.begin(); i != Empties.end(); ++i)
    {
	cout << "Empty" << i - Empties.begin() << ":\n";
	cout << *i;
    }
    cout << endl;

    cout << "===== Hexagons =====" << endl;
    for (Population::iterator i = Hexagons.begin(); i != Hexagons.end(); ++i)
    {
	cout << "Hexagon" << i - Hexagons.begin() << ":\n";
	cout << *i;
    }
    cout << endl;
    return EXIT_SUCCESS;
}
