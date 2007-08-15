/***********************************************************************
* Short Title: test of BGAlib.cpp
*
* Comments: test of Molecule::Rotate()
*
* $Id$
***********************************************************************/

#include <unistd.h>
#include "ParseArgs.hpp"
#include "Exceptions.hpp"
#include "BGAlib.hpp"

const int EXIT_INPUT_ERROR = 2;

using namespace std;

struct RunParameters
{
    RunParameters(int argc, char* argv[]);
    Molecule* pmol;
    // IO parameters
    string distfile;
    string inistru;
    string outstru;
    // Relaxation parameters
    double tol_dd;
private:
    void process_args(int argc, char* argv[]);
    void print_pars(ParseArgs& a);
    void print_help(ParseArgs& a);
    string version_string(string quote = "");
    list<string> validpars;
};

RunParameters::RunParameters(int argc, char* argv[])
{
    validpars.push_back("outstru");
    validpars.push_back("tol_dd");
    process_args(argc, argv);
}

void RunParameters::print_help(ParseArgs& a)
{
    // /usage:/;/;/-s/.*/"&\\n"/
    // /cou/;/;/s/^\s*"\(.*\)\\n"/\1/ | '[put! ='/*' | /;/put ='*/'
    cout << 
"usage: " << a.cmd_t << " distfile inistru.xyz [par1=val1 par2=val2...]\n"
"Options:\n"
"  -h, --help            display this message\n"
"  -V, --version         show program version\n"
"IO parameters:\n"
"  outstru=FILE          where to save relaxed molecule\n"
"Relaxation parameters:\n"
"  tol_dd=double         [Inf] distance is not used when dd=|d-d0|>=tol_dd\n"
"                        tol_dd=0.0 sets infinite distance table\n"
;
}

string RunParameters::version_string(string quote)
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

void RunParameters::print_pars(ParseArgs& a)
{
    // print out all run parameters
    // intro messages
    string hashsep(72, '#');
    cout << hashsep << endl;
    cout << "# " << a.cmd_t << endl;
    cout << version_string("# ");
    char hostname[255];
    gethostname(hostname, 255);
    cout << "# " << hostname << endl;
    time_t cur_time = time(NULL);
    cout << "# " << ctime(&cur_time);
    cout << hashsep << endl;
    bool second_hashsep = false;
    // outstru
    if (a.ispar("outstru"))
    {
        cout << "outstru=" << outstru << endl;
	second_hashsep = true;
    }
    // tol_dd
    if (a.ispar("tol_dd"))
    {
	cout << "tol_dd=" << tol_dd << endl;
	second_hashsep = true;
    }
    if (second_hashsep)
    {
	cout << hashsep << endl << endl;
    }
}

void RunParameters::process_args(int argc, char* argv[])
{
    char *short_options =
        "hV";
    // parameters and options
    option long_options[] = {
        {"help", 0, 0, 'h'},
        {"version", 0, 0, 'V'},
        {0, 0, 0, 0}
    };
    ParseArgs a(argc, argv, short_options, long_options);
    try {
        a.Parse();
    }
    catch (ParseArgsError) {
        exit(EXIT_INPUT_ERROR);
    }
    if (a.isopt("h") || argc == 1)
    {
        print_help(a);
        exit(EXIT_SUCCESS);
    }
    else if (a.isopt("V"))
    {
	cout << version_string();
        exit(EXIT_SUCCESS);
    }
    try {
	a.ValidatePars(validpars);
    }
    catch (ParseArgsError(e)) {
	cerr << e.what() << endl;
	exit(EXIT_INPUT_ERROR);
    }
    // assign run parameters
    // distfile, inistru
    switch (a.args.size())
    {
	case 0:
	    cerr << "distfile not specified" << endl;
	    exit(EXIT_INPUT_ERROR);
	case 1:
	    cerr << "initial structure not specified" << endl;
	    exit(EXIT_INPUT_ERROR);
	case 2:
	    distfile = a.args[0];
	    inistru = a.args[1];
	    break;
	default:
	    cerr << "too many command line arguments" << endl;
	    exit(EXIT_INPUT_ERROR);
    }
    // outstru
    if (a.ispar("outstru"))
    {
        outstru = a.pars["outstru"];
    }
    // relaxation parameters
    try {
	// tol_dd
	tol_dd = a.GetPar<double>("tol_dd", DOUBLE_MAX);
        DistanceTable dtab(distfile.c_str());
	pmol = new Molecule(dtab);
	pmol->tol_dd = tol_dd;
	if (tol_dd == 0.0)
	{
	    pmol->setMaxNAtoms(numeric_limits<int>().max());
	    pmol->ReadXYZ(inistru.c_str());
	    pmol->setMaxNAtoms(pmol->NAtoms());
	}
	else
	{
	    pmol->ReadXYZ(inistru.c_str());
	}
    }
    catch (IOError) {
        exit(EXIT_INPUT_ERROR);
    }
    catch (InvalidMolecule) {
        exit(EXIT_INPUT_ERROR);
    }
    catch (ParseArgsError(e)) {
	cerr << e.what() << endl;
	exit(EXIT_INPUT_ERROR);
    }
    // done
    print_pars(a);
}

list<int> worst_to_best_atom(Molecule& mol)
{
    typedef pair<double,int> AtomBadnessAndIdx;
    AtomBadnessAndIdx bi[mol.NAtoms()];
    for (int i = 0; i != mol.NAtoms(); ++i)
    {
	bi[i].first = mol.getAtom(i).Badness();
	bi[i].second = i;
    }
    sort(bi, bi + mol.NAtoms());
    // build list in reverse order
    list<int> w2b;
    for (AtomBadnessAndIdx* pbi = bi; pbi != bi + mol.NAtoms(); ++pbi)
    {
	w2b.push_front(pbi->second);
    }
    return w2b;
}

int main(int argc, char* argv[])
{
    RunParameters rp(argc, argv);
    Molecule& mol = *(rp.pmol);
    double mnb0 = DOUBLE_MAX;
    cout << "0 BC " << mol.NormBadness() << endl;
    for (int iteration = 1; mol.NormBadness() < mnb0; ++iteration)
    {
	list<int> indices = worst_to_best_atom(mol);
	for (   list<int>::iterator ii = indices.begin();
		ii != indices.end();  ++ii )
	{
	    mol.RelaxAtom(*ii);
	    cout << iteration << " R " << *ii
		<< " " << mol.NormBadness() << endl;
	}
	cout << iteration << " BC " << mol.NormBadness() << endl;
	if (eps_lt(mol.NormBadness(), mnb0))  mnb0 = mol.NormBadness();
    }
    if (!rp.outstru.empty())
    {
	try {
	    mol.WriteFile(rp.outstru.c_str());
	}
	catch (IOError(e)) {
	    cerr << e.what() << endl;
	    exit(EXIT_INPUT_ERROR);
	}
    }
    else
    {
	cout << endl;
	cout << "Relaxed molecule:\n" << endl;
	cout << mol;
    }
    return EXIT_SUCCESS;
}
