#include "GapsDispatcher.h"
#include "math/SIMD.h"

#include <Rcpp.h>

static RowMatrix convertRMatrix(const Rcpp::NumericMatrix &rmat)
{
    RowMatrix mat(rmat.nrow(), rmat.ncol());
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            mat(i,j) = rmat(i,j);
        }
    }
    return mat;
}

template <class Matrix>
static Rcpp::NumericMatrix createRMatrix(const Matrix &mat)
{
    Rcpp::NumericMatrix rmat(mat.nRow(), mat.nCol());
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            rmat(i,j) = mat(i,j);
        }
    }
    return rmat;
}

template <class T>
static Rcpp::List cogapsRun(T data, unsigned nPatterns,
unsigned maxIter, unsigned outputFrequency, unsigned seed, float alphaA,
float alphaP, float maxGibbsMassA, float maxGibbsMassP, bool messages,
bool singleCell, unsigned nCores)
{
    GapsDispatcher dispatcher(seed);

    dispatcher.setNumPatterns(nPatterns);
    dispatcher.setMaxIterations(maxIter);
    dispatcher.setOutputFrequency(outputFrequency);
    
    dispatcher.setAlpha(alphaA, alphaP);
    dispatcher.setMaxGibbsMass(maxGibbsMassA, maxGibbsMassP);

    dispatcher.printMessages(messages);
    dispatcher.singleCell(singleCell);
    dispatcher.setNumCoresPerSet(nCores);
    
    dispatcher.loadData(data);

    GapsResult result(dispatcher.run());
    return Rcpp::List::create(
        Rcpp::Named("Amean") = createRMatrix(result.Amean),
        Rcpp::Named("Pmean") = createRMatrix(result.Pmean),
        Rcpp::Named("Asd") = createRMatrix(result.Asd),
        Rcpp::Named("Psd") = createRMatrix(result.Psd)
    );
}

// [[Rcpp::export]]
Rcpp::List cogapsFromFile_cpp(const std::string &D, unsigned nPatterns,
unsigned maxIter, unsigned outputFrequency, unsigned seed, float alphaA,
float alphaP, float maxGibbsMassA, float maxGibbsMassP, bool messages,
bool singleCellRNASeq)
{
    return cogapsRun(D, nPatterns, maxIter, outputFrequency, seed,
        alphaA, alphaP, maxGibbsMassA, maxGibbsMassP, messages,
        singleCellRNASeq, 1);
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &D, unsigned nPatterns,
unsigned maxIter, unsigned outputFrequency, unsigned seed, float alphaA,
float alphaP, float maxGibbsMassA, float maxGibbsMassP, bool messages,
bool singleCellRNASeq, unsigned nCores)
{
    return cogapsRun(convertRMatrix(D), nPatterns, maxIter, outputFrequency, seed,
        alphaA, alphaP, maxGibbsMassA, maxGibbsMassP, messages,
        singleCellRNASeq, nCores);
}

/*
Rcpp::List cogapsFromCheckpoint_cpp(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, unsigned nEquil,
unsigned nSample, const std::string &cptFile)
{   
    GapsRunner runner(D, S, nFactor, nEquil, nSample, cptFile);
    return runner.run();
}
*/

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

#ifdef __GAPS_OPENMP__
    std::string openmp = "Compiled with OpenMP\n";
#else
    std::string openmp = "Compiler did not support OpenMP\n";
#endif

    return compiler + simd + openmp;
}
