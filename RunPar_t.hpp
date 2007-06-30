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
#include <valarray>
#include "BGAlib.hpp"
#include "ParseArgs.hpp"
#include "RegisterSVNId.hpp"

namespace {
RegisterSVNId RunPar_t_hpp_id("$Id$");
}

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

namespace VerboseFlag
{
    enum VerboseFlag { AD, WC, BC, AV, TS, ALL, VERBOSE_SIZE };
}

struct RunPar_t
{
    RunPar_t();
    void processArguments(int argc, char * argv[]);
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
    std::valarray<bool> verbose;
    bool verbose_mute;
    // Liga parameters
    size_t ndim;
    bool crystal;
    std::vector<double> latpar;
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

    std::string version_string(std::string quote="");
    std::list<std::string> validpars;
    std::vector<std::string> verbose_flag;

    void print_help(ParseArgs& a);
    void print_pars(ParseArgs& a);
    void fill_validpars();
    void fill_verbose();

};

#endif	// RUNPAR_T_HPP_INCLUDED
