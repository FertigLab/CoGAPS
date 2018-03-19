#include "Archive.h"
#include "GapsRunner.h"
#include "Random.h"
#include "SIMD.h"

#include <Rcpp.h>

// used to convert defined macro values into strings
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, unsigned nEquil,
unsigned nEquilCool, unsigned nSample, unsigned nOutputs, unsigned nSnapshots,
float alphaA, float alphaP, float maxGibbmassA, float maxGibbmassP, int seed,
bool messages, bool singleCellRNASeq, char whichMatrixFixed,
const Rcpp::NumericMatrix &FP, unsigned checkpointInterval,
const std::string &cptFile, unsigned pumpThreshold, unsigned nPumpSamples)
{
    // get seed, TODO do this on R side, multiple benefits (same seed in R, C++)
    uint32_t seedUsed = static_cast<uint32_t>(seed);
    if (seed < 0)
    {
        bpt::ptime epoch(boost::gregorian::date(1970,1,1));
        bpt::time_duration diff = bpt_now() - epoch;
        seedUsed = static_cast<uint32_t>(diff.total_milliseconds() % 1000);
    }

    // create internal state from parameters and run from there
    GapsRunner runner(D, S, nFactor, nEquil, nEquilCool, nSample,
        nOutputs, nSnapshots, alphaA, alphaP, maxGibbmassA, maxGibbmassP, seed,
        messages, singleCellRNASeq,  checkpointInterval, cptFile,
        whichMatrixFixed, FP);
    return runner.run();
}

// TODO add checksum to verify D,S matrices
// [[Rcpp::export]]
Rcpp::List cogapsFromCheckpoint_cpp(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, unsigned nEquil,
unsigned nSample, const std::string &fileName, const std::string &cptFile)
{   
    GapsRunner runner(D, S, nFactor, nEquil, nSample, cptFile);
    return runner.run();
}

// [[Rcpp::export]]
void displayBuildReport_cpp()
{
#if defined( __clang__ )
    Rcpp::Rcout << "Compiled with Clang\n";
#elif defined( __INTEL_COMPILER )
    Rcpp::Rcout << "Compiled with Intel ICC/ICPC\n";
#elif defined( __GNUC__ )
    Rcpp::Rcout << "Compiled with GCC v" << STR( __GNUC__ ) << "."
        << STR( __GNUC_MINOR__ ) << '\n';
#elif defined( _MSC_VER )
    Rcpp::Rcout << "Compiled with Microsoft Visual Studio\n";
#endif

#if defined( __GAPS_AVX__ )
    Rcpp::Rcout << "AVX enabled\n";
#elif defined( __GAPS_SSE__ )
    Rcpp::Rcout << "SSE enabled\n";
#else
    Rcpp::Rcout << "SIMD not enabled\n";
#endif
}
