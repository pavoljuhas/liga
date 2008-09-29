#include "EmbedPython.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////


void initializePython()
{
    static bool is_initialized = false;
    if (is_initialized)
    {
        return;
    }
    static int py_argc = 1;
    static char arg0[10] = "python";
    static char* py_argv[] = {arg0};
    Py_Initialize();
    PySys_SetArgv(py_argc, py_argv);
    is_initialized = true;
}


// End of file
