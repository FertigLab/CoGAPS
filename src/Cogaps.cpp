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

static bool isNull(const RowMatrix &mat)
{
    return mat.nRow() == 1 && mat.nCol() == 1;
}

static bool isNull(const std::string &path)
{
    return path.empty();
}

template <class T>
static Rcpp::List cogapsRun(const T &data, Rcpp::S4 params, const T &unc,
const RowMatrix &fixedMatrix, const std::string &checkpointInFile)
{
    GapsDispatcher dispatcher;

    // check if we're initializing with a checkpoint or not
    if (!isNull(checkpointInFile))
    {
        dispatcher.initialize(data, checkpointInFile);
    }
    else
    {
        dispatcher.initialize(data, params.slot("nPatterns"), params.slot("seed"));

        // set optional parameters
        dispatcher.setMaxIterations(params.slot("nIterations"));
        dispatcher.setOutputFrequency(params.slot("outputFrequency"));

        dispatcher.setSparsity(params.slot("alphaA"), params.slot("alphaP"),
            params.slot("singleCell"));
        
        dispatcher.setMaxGibbsMass(params.slot("maxGibbsMassA"),
            params.slot("maxGibbsMassP"));

        dispatcher.printMessages(params.slot("messages"));

        // check if running with a fixed matrix
        if (!isNull(fixedMatrix))
        {
            dispatcher.setFixedMatrix(params.slot("whichMatrixFixed"), fixedMatrix);
        }
    }

    // set the uncertainty matrix
    if (!isNull(unc))
    {
        dispatcher.setUncertainty(unc);
    }
    
    // set parameters that aren't saved in the checkpoint
    dispatcher.setNumCoresPerSet(params.slot("nCores"));
    dispatcher.setCheckpointInterval(params.slot("checkpointInterval"));
    dispatcher.setCheckpointOutFile(params.slot("checkpointOutFile"));

    // run the dispatcher and return the GapsResult in an R list
    GapsResult result(dispatcher.run());
    GAPS_ASSERT(result.meanChiSq > 0.f);
    return Rcpp::List::create(
        Rcpp::Named("Amean") = createRMatrix(result.Amean),
        Rcpp::Named("Pmean") = createRMatrix(result.Pmean),
        Rcpp::Named("Asd") = createRMatrix(result.Asd),
        Rcpp::Named("Psd") = createRMatrix(result.Psd),
        Rcpp::Named("seed") = result.seed,
        Rcpp::Named("meanChiSq") = result.meanChiSq,
        Rcpp::Named("diagnostics") = Rcpp::List::create()
    );
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp_from_file(const std::string &data,
const Rcpp::List &allParams,
const std::string &uncertainty,
const Rcpp::NumericVector &indices=Rcpp::NumericVector(),
const Rcpp::NumericMatrix &initialA=Rcpp::NumericMatrix(),
const Rcpp::NumericMatrix &initialP=Rcpp::NumericMatrix())
{
    return cogapsRun(data, params, unc, convertRMatrix(fixedMatrix),
        checkpointInFile);
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &data,
const Rcpp::List &allParams,
const Rcpp::NumericMatrix &uncertainty,
const Rcpp::NumericVector &indices=Rcpp::NumericVector(),
const Rcpp::NumericMatrix &initialA=Rcpp::NumericMatrix(),
const Rcpp::NumericMatrix &initialP=Rcpp::NumericMatrix())
{
    return cogapsRun(convertRMatrix(data), params, convertRMatrix(unc),
        convertRMatrix(fixedMatrix), checkpointInFile);
}

// [[Rcpp::export]]
std::string getBuildReport_cpp()
{
    return buildReport();
}
