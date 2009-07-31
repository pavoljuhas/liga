/***********************************************************************
* Short Title: run parameters for gizaliga application
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef RUNPAR_T_HPP_INCLUDED
#define RUNPAR_T_HPP_INCLUDED

#include <deque>
#include <string>
#include <memory>
#include "ParseArgs.hpp"
#include "TraceId_t.hpp"
#include "Molecule.hpp"

struct RunPar_t
{
    RunPar_t(int argc, char* argv[]);
    void processArguments(int argc, char* argv[]);
    // parsed input arguments
    std::auto_ptr<ParseArgs> args;
    // Output option
    bool trace;
    // IO parameters
    std::string distfile;
    std::string inistru;
    std::string outstru;
    std::string outfmt;
    int saverate;
    bool saveall;
    std::string frames;
    int framesrate;
    std::deque<TraceId_t> framestrace;
    std::vector<bool> verbose;
    // Liga parameters
    size_t ndim;
    bool crystal;
    std::vector<double> latpar;
    double rmax;
    bool distreuse;
    double tolcost;
    ChemicalFormula formula;
    std::vector<int> fixed_atoms;
    double maxcputime;
    int rngseed;
    double promotefrac;
    bool promoterelax;
    bool demoterelax;
    int ligasize;
    double stopgame;
    int seasontrials;
    std::string trialsharing;
    // generated data
    std::auto_ptr<Molecule> mol;
    int base_level;
    // Constrains
    std::vector<double> bangle_range;
    double max_dist;

private:

    std::string version_string(std::string quote="");
    const std::list<std::string>& validpars() const;
    const std::string& joined_verbose_flags() const;

    void print_help();
    void print_pars();
    void fill_validpars();

};

#endif	// RUNPAR_T_HPP_INCLUDED
