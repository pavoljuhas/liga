/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test of push_good_pyramids()
*
* $Id$
***********************************************************************/

#include "BGAlib.hpp"

int main()
{
    using namespace std;

    // line distance table
    valarray<double> pdT(sqrt(2.0), 6);
    SandSphere ssT(100, 6, &pdT[0]);
    ssT.SetGridTol(1.1);
    // tetrahedron coordinates
    double px0[4] = {0, 1, 1, 0};
    double py0[4] = {0, 1, 0, 1};
    double pz0[4] = {0, 0, 1, 1};
    Molecule tetrahedron(&ssT, 4, px0, py0, pz0);
    // evaluate and print fitness for line
    cout << "tetrahedron:\n";
    cout << tetrahedron;
    tetrahedron.PrintBadness();
    tetrahedron.PrintFitness();
    cout << endl;

    cout << "tetrahedron.Pop(0):" << endl;
    tetrahedron.Pop(0);
    cout << tetrahedron;
    tetrahedron.PrintBadness();
    tetrahedron.PrintFitness();
    cout << endl;

    cout << "trying to build tetrahedron:" << endl;
    tetrahedron.Clear();
    tetrahedron.Add(0,0,0);
    cout << tetrahedron;
    tetrahedron.PrintBadness();

    cout << "push_good_distances" << endl;
    double afit[4] = {1,1,1,1};
    vector<Atom_t> v1;
    tetrahedron.push_good_distances(v1, afit, 1);
    tetrahedron.Add(v1[0]);
    cout << tetrahedron;
    tetrahedron.PrintBadness();
    tetrahedron.PrintFitness();

    cout << "push_good_triangles" << endl;
    vector<Atom_t> v2;
    tetrahedron.push_good_triangles(v2, afit, 1);
    tetrahedron.Add(v2[0]);
    cout << tetrahedron;
    tetrahedron.PrintBadness();
    tetrahedron.PrintFitness();

    cout << "push_good_pyramids" << endl;
    vector<Atom_t> v3;
    tetrahedron.push_good_pyramids(v3, afit, 1);
    tetrahedron.Add(v3[0]);
    cout << tetrahedron;
    tetrahedron.PrintBadness();
    tetrahedron.PrintFitness();

    /*
    tetrahedron.push_good_pyramids(vertices, afit, 8);
    for (int i=0; i<vertices.size(); ++i)
    {
	cout << 'V' << i << ": [ " <<
	    vertices[i].h << ' ' <<
	    vertices[i].k << ' ' <<
	    vertices[i].l << ' ' <<
	    "]  " << vertices[i].Badness() << endl;
    }
    */

    return EXIT_SUCCESS;
}
