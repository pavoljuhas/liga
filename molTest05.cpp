/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test Molecule::ABadnessAt()
*     molTest05 m writes matlab code for plotting badness matrix
*
* $Id$
***********************************************************************/

#include <unistd.h>
#include "BGAlib.hpp"

int main(int argc, char *argv[])
{
    using namespace std;

    // distance table
    SandSphere ssHex(100, "hexagon.dss");

    Molecule mol1(&ssHex);
    mol1.ReadXY("hexagon.dxy");

    bool showMatlab = 
	(argc > 1 && strlen(argv[1]) > 0 && tolower(argv[1][0]) == 'm');
    if (!showMatlab)
    {
	cout << "mol1:" << endl;
	cout << mol1; mol1.PrintBadness(); cout << endl;
	cout << "mol1.Pop(0):" << endl;
	mol1.Pop(0);
	cout << mol1; mol1.PrintBadness(); cout << endl;
	cout << "mol1.ABadnessAt(0, 87) = " << mol1.ABadnessAt(0, 87) << endl;
	cout << "mol1.ABadnessAt(50, 0) = " << mol1.ABadnessAt(50, 0) << endl;
	cout << "mol1.ABadnessAt(-50, 0) = " << mol1.ABadnessAt(-50, 0) << endl;
	cout << endl;
    }
    else
    {
	mol1.Pop(0);
	ostringstream oss;
	oss << "hk = [" << endl << mol1 << "];" << endl;
	string head = oss.str();
	for (   string::size_type sp = head.find('#', 0);
		sp != string::npos; sp = head.find('#', ++sp)
	    )
	{
	    head[sp] = '%';
	}
	cout << head;
	cout << "X = " << -ssHex.gridmax << ":" << ssHex.gridmax << "; " <<
	    "Y = X;" << endl;
	cout << "Z = [" << endl;
	for (int i = -ssHex.gridmax; i <= ssHex.gridmax; ++i)
	{
	    for (int j = -ssHex.gridmax; j <= ssHex.gridmax; ++j)
	    {
		cout << ' ' << mol1.ABadnessAt(i, j);
	    }
	    cout << endl;
	}
	cout << "]';" << endl;
	cout << "echo on" << endl;
	cout << "clf; imagesc(X,Y,Z,[0 5])" << endl;
	cout << "hold on; hs=plot(hk(:,1),hk(:,2),'r*');" << endl;
	cout << "echo off" << endl;
    }
    return EXIT_SUCCESS;
}
