/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: create SandSphere object using 3 different constructors
*
* $Id$
***********************************************************************/

#include "BGAlib.hpp"

int main()
{
    double pd[6] = {1.0, 1.0, 1.4142, 1.4142, 1.0, 1.0};
    vector<double> vd(6);
    vd[0] = 1.0;
    vd[1] = 1.0;
    vd[2] = 1.4142;
    vd[3] = 1.4142;
    vd[4] = 1.0;
    vd[5] = 1.0;
    char *distfile = "square.dss";
    SandSphere *ss1, *ss2, *ss3;
    try {
	ss1 = new SandSphere(100, vd);
	ss2 = new SandSphere(100, 6, pd);
	ss3 = new SandSphere(100, distfile);
    }
    catch (...) {
	return EXIT_FAILURE;
    }
    delete ss1, ss2, ss3;
    return EXIT_SUCCESS;
}
