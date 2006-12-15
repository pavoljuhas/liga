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
    // process arguments
    RunPar_t rp;
    try	{
	rp.processArguments(argc, argv);
    }
    catch (IOError(e)) {
	cerr << e.what() << endl;
	exit(EXIT_INPUT_ERROR);
    }
    catch (ParseArgsError(e)) {
	cerr << e.what() << endl;
	exit(EXIT_INPUT_ERROR);
    }
    Liga_t liga(&rp);
    // watch for HUP
    signal(SIGHUP, SIGHUP_handler);
    liga.stopFlag(&SIGHUP_received);
    // main loop
    for (liga.prepare(); !liga.finished(); liga.playSeason()) { }
    // figure out exit code
    int exit_code;
    if (SIGHUP_received)	    exit_code = SIGHUP + 128;
    else if (liga.solutionFound())  exit_code = EXIT_SUCCESS;
    else			    exit_code = EXIT_FAILURE;
    return exit_code;
}

// End of file
