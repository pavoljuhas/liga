/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test Molecule::Evolve(), Molecule::Degenerate()
*     molTest06 distance_file.dss [snapshot_file.dxy]
*
* $Id$
***********************************************************************/

#include <unistd.h>
#include "BGAlib.hpp"

int main(int argc, char *argv[])
{
    using namespace std;

    if (argc == 1)
    {
	cerr << "usage: " <<
	    "molTest06 distance_file.dss [snapshot_file.dxy]" << endl;
	return EXIT_SUCCESS;
    }
    // here argc > 1
    char *distance_file = argv[1];
    SandSphere* ss;
    try
    {
	ss = new SandSphere(100, distance_file);
    }
    catch (IOError)
    {
	return EXIT_FAILURE;
    }
    char *snapshot_file = NULL;
    if (argc > 2 && strlen(argv[2]) > 0)
    {
	snapshot_file = argv[2];
	cout << "molecule snapshots go to " << snapshot_file << endl;
    }
	
    Molecule mol(ss);

    cout << "use (e) to evolve, (d) to degenerate, (q) to quit\n\n";
    cout << "mol:" << endl << mol; mol.PrintBadness();
    if (snapshot_file != NULL)
    {
	mol.WriteXY(snapshot_file);
    }
    char c = '0';
    while (c != 'q')
    {
	c = tolower(cin.get());
	bool mol_changed = false;
	switch (c)
	{
	    case 'e':
		if (mol.NAtoms == mol.MaxAtoms)
		    cout << "full size molecule, Evolve() ignored...\n";
		else
		{
		    mol.Evolve();
		    mol_changed = true;
		}
		break;
	    case 'd':
		if (mol.NAtoms == 0)
		    cout << "empty size molecule, Degenerate() ignored...\n";
		else
		{
		    mol.Degenerate();
		    mol_changed = true;
		}
		break;
	}
	if (mol_changed)
	{
	    cout << "mol:" << endl << mol; mol.PrintBadness();
	    if (snapshot_file != NULL)
	    {
		mol.WriteXY(snapshot_file);
	    }
	}
    }
    return EXIT_SUCCESS;
}
