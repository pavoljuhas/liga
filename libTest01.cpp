/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test of random_choose_few()
*     libTest01 [repeats=1]
*
* $Id$
***********************************************************************/

#include "BGAlib.hpp"

int main(int argc, char* argv[])
{
    using namespace std;
    const int N = 10;
    int repeats = 1;
    if (argc > 1)
	repeats = atoi(argv[1]);
    for (int c = 0; c < repeats; ++c)
    {
	for (int i = 0; i <= N; ++i)
	{
	    vector<int> vec = random_choose_few(i, N);
	    cout << i << " of " << N << ':';
	    for (vector<int>::iterator vi = vec.begin(); vi != vec.end(); ++vi)
		cout << ' ' << *vi;
	    cout << endl;
	}
    }
    return EXIT_SUCCESS;
}
