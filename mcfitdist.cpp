/***********************************************************************
* Short Title: MC molecule reconstruction from distance table
*
* Comments: Metropolis algorithm
*
* $Id$
***********************************************************************/

#include <limits>
#include <unistd.h>
#include <signal.h>
#include "ParseArgs.hpp"
#include "BGAlib.hpp"

typedef valarray<Atom_t> VAA;

////////////////////////////////////////////////////////////////////////
// RunPar_t
////////////////////////////////////////////////////////////////////////

struct RunPar_t
{
    RunPar_t();
    VAA ProcessArguments(int argc, char * argv[]);
    // IO parameters
    string distfile;
    string inistru;
    string outstru;
    int saverate;
    string frames;
    int framesrate;
    // MC parameters
    double boxsize;
    double kbt;
    double tol_bad;
    int seed;
    string penalty;
private:
    void print_help(ParseArgs& a);
    string version_string(string quote = "");
    list<string> validpars;
};

RunPar_t::RunPar_t()
{
    char *pnames[] = {
	"distfile", "inistru",
	"outstru", "saverate", "frames", "framesrate",
	"boxsize", "kbt", "tol_bad", "seed",
	"penalty" };
    validpars.insert(validpars.end(),
	    pnames, pnames+sizeof(pnames)/sizeof(char*));
}

void RunPar_t::print_help(ParseArgs& a)
{
    // /usage:/;/;/-s/.*/"&\\n"/
    // /cou/;/;/s/^\s*"\(.*\)\\n"/\1/ | '[put! ='/*' | /;/put ='*/'
    cout << 
"usage: " << a.cmd_t << "[-p PAR_FILE] [DISTFILE] [par1=val1 par2=val2...]\n"
"run MC molecule simulation using distances from DISTFILE.  Parameters can\n"
"be set in PAR_FILE or on the command line, which overrides PAR_FILE.\n"
"Options:\n"
"  -p, --parfile=FILE    read parameters from FILE\n"
"  -h, --help            display this message\n"
"  -v, --version         show program version\n"
"IO parameters:\n"
"  distfile=FILE         target distance table\n"
"  inistru=FILE          initial structure [empty box]\n"
"  outstru=FILE          where to save the best full molecule\n"
"  saverate=int          [10] minimum iterations between outstru updates\n"
"  frames=FILE           save intermediate structures to FILE.liga_round\n"
"  framesrate=int        [10] number of iterations between frame saves\n"
"Liga parameters\n"
"  boxsize=double        [0.1] size of box of possible MC step\n"
"  kbt=double            [0.001] Boltzman factor\n"
"  tol_bad=double        [1E-4] target normalized molecule badness\n"
"  seed=int              seed of random number generator\n"
"  penalty=string        dd penalty function [pow2], fabs, well\n"
;
}

string RunPar_t::version_string(string quote)
{
    using namespace std;
    ostringstream oss;
    oss << quote
        << "$Id$" << endl
#   if defined(__DATE__) && defined(__TIME__)
	<< quote << "compiled " __DATE__ " " __TIME__ << endl
#   endif
        ;
    return oss.str();
}

