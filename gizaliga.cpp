/***********************************************************************
* Short Title: molecule reconstruction from distance table
*
* Comments: self-tuning competition search for molecular configuration
*
* $Id$
***********************************************************************/

#include <limits>
#include <unistd.h>
#include <signal.h>
#include "ParseArgs.hpp"
#include "BGAlib.hpp"
#include "Liga_t.hpp"

const int EXIT_INPUT_ERROR = 2;

////////////////////////////////////////////////////////////////////////
// SIGHUP handling
////////////////////////////////////////////////////////////////////////

int SIGHUP_received = 0;

void SIGHUP_handler(int signum)
{
    // die on 2nd call
    if (SIGHUP_received)    exit(128 + signum);
    SIGHUP_received = signum;
}

////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    RunPar_t* rp = NULL;
    Liga_t* liga = NULL;
    // Catch exceptions
    try	{
	// process arguments
	rp = new RunPar_t();
	rp->processArguments(argc, argv);
	liga = new Liga_t(rp);
	// watch for HUP
	signal(SIGHUP, SIGHUP_handler);
	liga->useStopFlag(&SIGHUP_received);
	// main loop
	liga->prepare();
	while (!liga->finished())    liga->playSeason();
	liga->printSummary();
    }
    catch (IOError(e)) {
	cerr << e.what() << endl;
	delete liga; delete rp;
	exit(EXIT_INPUT_ERROR);
    }
    catch (ParseArgsError(e)) {
	cerr << e.what() << endl;
	delete liga; delete rp;
	exit(EXIT_INPUT_ERROR);
    }
    catch (runtime_error(e)) {
	cerr << e.what() << endl;
	delete liga; delete rp;
	exit(EXIT_INPUT_ERROR);
    }
    // figure out exit code
    int exit_code;
    if (SIGHUP_received)	    exit_code = SIGHUP + 128;
    else if (liga->solutionFound()) exit_code = EXIT_SUCCESS;
    else			    exit_code = EXIT_FAILURE;
    // clean up
    delete liga; delete rp;
    return exit_code;
}

// End of file
