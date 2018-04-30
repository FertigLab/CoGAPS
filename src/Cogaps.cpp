#include "SIMD.h"
#include "GapsRunner.h"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, unsigned nEquil,
unsigned nEquilCool, unsigned nSample, unsigned nOutputs, unsigned nSnapshots,
float alphaA, float alphaP, float maxGibbmassA, float maxGibbmassP,
unsigned seed, bool messages, bool singleCellRNASeq, char whichMatrixFixed,
const Rcpp::NumericMatrix &FP, unsigned checkpointInterval,
const std::string &cptFile, unsigned pumpThreshold, unsigned nPumpSamples)
{
    // create internal state from parameters and run from there
    GapsRunner runner(D, S, nFactor, nEquil, nEquilCool, nSample,
        nOutputs, nSnapshots, alphaA, alphaP, maxGibbmassA, maxGibbmassP, seed,
        messages, singleCellRNASeq, checkpointInterval, cptFile,
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

// used to convert defined macro values into strings
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

// [[Rcpp::export]]
std::string getBuildReport_cpp()
{
#if defined( __clang__ )
    std::string compiler = "Compiled with Clang\n";
#elif defined( __INTEL_COMPILER )
    std::string compiler = "Compiled with Intel ICC/ICPC\n";
#elif defined( __GNUC__ )
    std::string compiler = "Compiled with GCC v" + std::string(STR( __GNUC__ ))
    + "." + std::string(STR( __GNUC_MINOR__ )) + '\n';
#elif defined( _MSC_VER )
    std::string compiler = "Compiled with Microsoft Visual Studio\n";
#endif

#if defined( __GAPS_AVX__ )
    std::string simd = "AVX enabled\n";
#elif defined( __GAPS_SSE__ )
    std::string simd = "SSE enabled\n";
#else
    std::string simd = "SIMD not enabled\n";
#endif
    return compiler + simd;
}
