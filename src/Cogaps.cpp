#include "GapsDispatcher.h"
#include "math/SIMD.h"

#include <Rcpp.h>

#include <string>

static std::vector<unsigned> convertRVec(const Rcpp::NumericVector &rvec)
{
    std::vector<unsigned> vec;
    for (unsigned i = 0; i < rvec.size(); ++i)
    {
        vec.push_back(rvec[i]);
    }
    return vec;
}

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

static bool isNull(const Rcpp::NumericMatrix &mat)
{
    return mat.nrow() == 1 && mat.ncol();
}

static bool isNull(const Rcpp::NumericVector &vec)
{
    return vec.size() == 1;
}

static bool isNull(const std::string &path)
{
    return path.empty();
}

template <class T>
static Rcpp::List cogapsRun(const T &data, const Rcpp::List &allParams,
const T &uncertainty, const Rcpp::NumericVector &indices,
const Rcpp::NumericMatrix &initialA, const Rcpp::NumericMatrix &initialP)
{
    GapsDispatcher dispatcher;

    // check if we're initializing with a checkpoint or not
    if (!isNull(allParams["checkpointInFile"]))
    {
        dispatcher.initialize(data, allParams["transposeData"],
            allParams["checkpointInFile"]);
    }
    else
    {
        // check to see if we're subsetting the data
        if (!isNull(allParams["gaps"].slot("distributed")))
        {
            dispatcher.initialize(data, allParams["transposeData"],
                allParams["gaps"].slot("distributed") == "genome-wide",
                convertRVec(indices), allParams["gaps"].slot("nPatterns"),
                allParams["gaps"].slot("seed"));
        }
        else
        {
            dispatcher.initialize(data, allParams["transposeData"],
                allParams["gaps"].slot("nPatterns"), allParams["gaps"].slot("seed"));
        }

        // set optional parameters
        dispatcher.setMaxIterations(allParams["gaps"].slot("nIterations"));
        dispatcher.setSparsity(allParams["gaps"].slot("alphaA"),
            allParams["gaps"].slot("alphaP"), allParams["gaps"].slot("singleCell"));
        dispatcher.setMaxGibbsMass(allParams["gaps"].slot("maxGibbsMassA"),
            allParams["gaps"].slot("maxGibbsMassP"));

        // set initial values for A and P matrix
        if (!isNull(initialA))
        {
            dispatcher.setAMatrix(convertRMatrix(initialA));
        }
        if (!isNull(initialP))
        {
            dispatcher.setPMatrix(convertRMatrix(initialP));
        }

        // check if running with a fixed matrix
        if (!isNull(allParams["whichMatrixFixed"]))
        {
            dispatcher.setFixedMatrix(allParams["whichMatrixFixed"]);
        }
    }

    // set the uncertainty matrix
    if (!isNull(uncertainty))
    {
        dispatcher.setUncertainty(uncertainty);
    }
    
    // set parameters that aren't saved in the checkpoint
    dispatcher.setNumCoresPerSet(params.slot("nCores"));
    dispatcher.printMessages(params.slot("messages"));
    dispatcher.setOutputFrequency(params.slot("outputFrequency"));
    dispatcher.setCheckpointOutFile(params.slot("checkpointOutFile"));
    dispatcher.setCheckpointInterval(params.slot("checkpointInterval"));

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
const Rcpp::NumericVector &indices=Rcpp::NumericVector(1),
const Rcpp::NumericMatrix &initialA=Rcpp::NumericMatrix(1,1),
const Rcpp::NumericMatrix &initialP=Rcpp::NumericMatrix(1,1))
{
    return cogapsRun(data, allParams, uncertainty, indices, initialA, initialP);
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &data,
const Rcpp::List &allParams,
const Rcpp::NumericMatrix &uncertainty,
const Rcpp::NumericVector &indices=Rcpp::NumericVector(1),
const Rcpp::NumericMatrix &initialA=Rcpp::NumericMatrix(1,1),
const Rcpp::NumericMatrix &initialP=Rcpp::NumericMatrix(1,1))
{
    return cogapsRun(convertRMatrix(data), allParams,
        convertRMatrix(uncertainty), indices, initialA, initialP);
}

// [[Rcpp::export]]
std::string getBuildReport_cpp()
{
    return buildReport();
}
