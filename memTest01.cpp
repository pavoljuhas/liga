/***********************************************************************
* Short Title: test liga memory usage
*
* Comments: just build a completely filled liga for a specified
*     DistanceTable and ligasize, and wait for <Enter> to quit
*
* $Id$
***********************************************************************/

#include <limits>
#include <unistd.h>
#include <signal.h>
#include "ParseArgs.hpp"
#include "BGAlib.hpp"

////////////////////////////////////////////////////////////////////////
// Division_t
////////////////////////////////////////////////////////////////////////
typedef Molecule* PMOL;
struct Division_t : public vector<Molecule*>
{
public:
    // constructors
    Division_t(int s)  : vector<PMOL>(), max_size(s) { }
    Division_t(const Division_t& div0) :
	vector<PMOL>(div0), max_size(div0.max_size) { }
    ~Division_t()
    {
	for (iterator ii = begin(); ii != end(); ++ii)
	    delete *ii;
    }
    Division_t& operator= (const vector<PMOL>& div0)
    {
	*this = div0;
	return *this;
    }
    Division_t& operator= (const Division_t& div0)
    {
	*this = vector<PMOL>(div0);
	max_size = div0.max_size;
	return *this;
    }
    int find_winner();
    PMOL& winner();
    int find_looser();
    PMOL& looser();
    PMOL& best();
    inline bool full() { return !(size() < max_size); }
    inline int fullsize() { return max_size; }
private:
    mutable int max_size;
};

int Division_t::find_winner()
{
    // evaluate molecule fitness
    valarray<double> vmfit(size());
    double *pd = &vmfit[0];
    for (iterator mi = begin(); mi != end(); ++mi, ++pd)
	*pd = (*mi)->NormBadness();
    // then get the reciprocal value
    vmfit = vdrecipw0(vmfit);
    double *mfit = &vmfit[0];
    int idx = random_wt_choose(1, mfit, size()).front();
    return idx;
}

PMOL& Division_t::winner()
{
    return at(find_winner());
}

int Division_t::find_looser()
{
    // evaluate molecule fitness
    valarray<double> vmbad(size());
    double *pd = &vmbad[0];
    for (iterator mi = begin(); mi != end(); ++mi, ++pd)
	*pd = (*mi)->NormBadness();
    double *mbad = &vmbad[0];
    int idx = random_wt_choose(1, mbad, size()).front();
    return idx;
}

PMOL& Division_t::looser()
{
    return at(find_looser());
}

bool comp_PMOL_Badness(const PMOL& lhs, const PMOL& rhs)
{
    return lhs->Badness() < rhs->Badness();
}

PMOL& Division_t::best()
{
    // evaluate molecule fitness
    iterator pm = min_element(begin(), end(), comp_PMOL_Badness);
    return *pm;
}


////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
	cout << "usage: memTest01 distfile ligasize" << endl;
	exit(EXIT_SUCCESS);
    }
    char* distfile = argv[1];
    int ligasize = atoi(argv[2]);
    DistanceTable dtab(distfile);
    Molecule mol(dtab);
    // initialize liga divisions, primitive divisions have only 1 team
    vector<Division_t> liga;
    for (int i = 0; i <= mol.max_NAtoms(); ++i)
    {
        int divsize = (i < 2) ? 1 : ligasize;
        liga.push_back(Division_t(divsize));
    }
    // put initial molecule to the zeroth division
    liga[0].push_back(new Molecule(mol));
    // fill all the remaining levels
    for (int level = 1; level <= mol.max_NAtoms(); ++level)
    {
	PMOL lower_team = liga[level-1].back();
	while(! liga[level].full())
	{
	    double xn = gsl_rng_uniform(BGA::rng) * dtab.back();
	    double yn = gsl_rng_uniform(BGA::rng) * dtab.back();
	    double zn = gsl_rng_uniform(BGA::rng) * dtab.back();
	    Molecule* new_team = new Molecule(*lower_team);
	    new_team->Add(xn, yn, zn);
	    liga[level].push_back(new_team);
	}
    }
    cout << "Filled liga system with ligasize = " << ligasize
	<< " natoms = " << mol.max_NAtoms() << endl
	<< "check memory usage with ps -C memTest01 -o rss,args" << endl
	<< "and then press <Enter> to exit" << endl;
    string s;
    getline(cin, s);
    exit(EXIT_SUCCESS);
}
