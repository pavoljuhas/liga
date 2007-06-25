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
#include "BGAlib.hpp"
#include "ParseArgs.hpp"

struct TraceId_t
{
    int season;
    int level;
    int index;
};

struct SeedClusterInfo
{
    int level;
    int number;
    int trials;
};

struct RunPar_t
{
    RunPar_t();
    void processArguments(int argc, char * argv[]);
    // Output option
    bool quiet;
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
    // Liga parameters
    size_t ndim;
    double tol_dd;
    double tol_bad;
    int natoms;
    std::vector<int> fixed_atoms;
    std::vector<SeedClusterInfo> seed_clusters;
    int centersize;
    double maxcputime;
    int rngseed;
    double evolve_frac;
    bool evolve_relax;
    bool degenerate_relax;
    int ligasize;
    double stopgame;
    int seasontrials;
    std::string trials_sharing;
    double lookout_prob;
    // generated data
    Molecule mol;
    int base_level;
    // Constrains
    std::vector<double> bangle_range;
    double max_dist;

private:

    void print_help(ParseArgs& a);
    std::string version_string(std::string quote = "");
    void print_pars(ParseArgs& a);
    std::list<std::string> validpars;
    void fill_validpars();

};

#endif	// RUNPAR_T_HPP_INCLUDED
