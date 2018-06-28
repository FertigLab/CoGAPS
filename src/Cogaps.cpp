#include "GapsDispatcher.h"
#include "math/SIMD.h"

#include <Rcpp.h>

#include <string>

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

static bool nonNullUncertainty(const RowMatrix &mat)
{
    return mat.nRow() > 1 || mat.nCol() > 1;
}

static bool nonNullUncertainty(const std::string &path)
{
    return !path.empty();
}

template <class T>
static Rcpp::List cogapsRun(const T &data, const T &unc, unsigned nPatterns,
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

    if (nonNullUncertainty(unc))
    {
        dispatcher.setUncertainty(unc);
    }

    GapsResult result(dispatcher.run());
    return Rcpp::List::create(
        Rcpp::Named("Amean") = createRMatrix(result.Amean),
        Rcpp::Named("Pmean") = createRMatrix(result.Pmean),
        Rcpp::Named("Asd") = createRMatrix(result.Asd),
        Rcpp::Named("Psd") = createRMatrix(result.Psd)
    );
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp_from_file(const std::string &data, const std::string &unc,
unsigned nPatterns, unsigned maxIterations, unsigned outputFrequency,
uint32_t seed, float alphaA, float alphaP, float maxGibbsMassA,
float maxGibbsMassP, bool messages, bool singleCell,
const std::string &checkpointOutFile, unsigned nCores)
{
    return cogapsRun(data, unc, nPatterns, maxIterations, outputFrequency, seed,
        alphaA, alphaP, maxGibbsMassA, maxGibbsMassP, messages, singleCell,
        nCores);
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &data,
const Rcpp::NumericMatrix &unc, unsigned nPatterns, unsigned maxIterations,
unsigned outputFrequency, uint32_t seed, float alphaA, float alphaP,
float maxGibbsMassA, float maxGibbsMassP, bool messages, bool singleCell,
const std::string &checkpointOutFile, unsigned nCores)
{
    return cogapsRun(convertRMatrix(data), convertRMatrix(unc), nPatterns,
        maxIterations, outputFrequency, seed, alphaA, alphaP, maxGibbsMassA,
        maxGibbsMassP, messages, singleCell, nCores);
}

// [[Rcpp::export]]
std::string getBuildReport_cpp()
{
    return buildReport();
}
