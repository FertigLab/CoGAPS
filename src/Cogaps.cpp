#include "GapsRunner.h"
#include "math/SIMD.h"

#include <Rcpp.h>

#include <string>

// these are helper functions for converting matrix/vector types
// to and from R objects

static std::vector<unsigned> convertRVec(const Rcpp::NumericVector &rvec)
{
    std::vector<unsigned> vec;
    for (unsigned i = 0; i < rvec.size(); ++i)
    {
        vec.push_back(rvec[i]);
    }
    return vec;
}

static Matrix convertRMatrix(const Rcpp::NumericMatrix &rmat)
{
    Matrix mat(rmat.nrow(), rmat.ncol());
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

// this provides a standard way for communicating which parameters
// are null between R and C++

static bool isNull(const Matrix &mat)
{
    return mat.nRow() == 1 && mat.nCol() == 1;
}

static bool isNull(const Rcpp::NumericMatrix &mat)
{
    return mat.nrow() == 1 && mat.ncol() == 1;
}

static bool isNull(const Rcpp::NumericVector &vec)
{
    return vec.size() == 1;
}

static bool isNull(const std::string &path)
{
    return path.empty();
}

// this is the main function that creates a GapsRunner and runs CoGAPS

template <class DataType>
static Rcpp::List cogapsRun(const DataType &data, const Rcpp::List &allParams,
const DataType &uncertainty, const Rcpp::NumericVector &indices,
const Rcpp::NumericMatrix &fixedMatrix)
{
    // convert string parameters
    Rcpp::S4 gapsParams = allParams["gaps"];
    std::string checkpointInFile = Rcpp::as<std::string>(allParams["checkpointInFile"]);
    std::string distributed = Rcpp::as<std::string>(gapsParams.slot("distributed"));
    GAPS_ASSERT(distributed == "genome-wide" || distributed == "single-cell");
    bool partitionRows = (distributed == "genome-wide");
    
    // read number of patterns from checkpoint file
    unsigned nPatterns = gapsParams.slot("nPatterns");
    unsigned seed = gapsParams.slot("seed"); // so we can return seed
    Archive ar(checkpointInFile, ARCHIVE_READ);
    if (!isNull(checkpointInFile))
    {
        gaps::random::load(ar);
        ar >> nPatterns >> seed;
    }

    // construct GapsRunner
    GapsRunner runner(data, allParams["transposeData"], nPatterns, seed,
        partitionRows, convertRVec(indices));

    // populate GapsRunner from checkpoint file
    if (!isNull(checkpointInFile))
    {
        ar >> runner;
        ar.close();
    }
    else
    {
        // set fixed matrix
        if (!isNull(fixedMatrix))
        {
            std::string which = Rcpp::as<std::string>(allParams["whichMatrixFixed"]);
            GAPS_ASSERT(!isNull(which));
            runner.setFixedMatrix(which[0], convertRMatrix(fixedMatrix));
        }

        // set parameters that would be saved in the checkpoint 
        gaps::random::setSeed(seed);
        runner.setMaxIterations(gapsParams.slot("nIterations"));
        runner.setSparsity(gapsParams.slot("alphaA"),
            gapsParams.slot("alphaP"), gapsParams.slot("singleCell"));
        runner.setMaxGibbsMass(gapsParams.slot("maxGibbsMassA"),
            gapsParams.slot("maxGibbsMassP"));
    }

    // set uncertainty
    if (!isNull(uncertainty))
    {
        runner.setUncertainty(uncertainty, allParams["transposeData"],
            partitionRows, convertRVec(indices));
    }
    
    // set parameters that aren't saved in the checkpoint
    runner.setMaxThreads(allParams["nThreads"]);
    runner.setPrintMessages(allParams["messages"]);
    runner.setOutputFrequency(allParams["outputFrequency"]);
    runner.setCheckpointOutFile(allParams["checkpointOutFile"]);
    runner.setCheckpointInterval(allParams["checkpointInterval"]);

    // run cogaps and return the GapsResult in an R list
    GapsResult result(runner.run());
    GAPS_ASSERT(result.meanChiSq > 0.f);
    return Rcpp::List::create(
        Rcpp::Named("Amean") = createRMatrix(result.Amean),
        Rcpp::Named("Pmean") = createRMatrix(result.Pmean),
        Rcpp::Named("Asd") = createRMatrix(result.Asd),
        Rcpp::Named("Psd") = createRMatrix(result.Psd),
        Rcpp::Named("seed") = seed,
        Rcpp::Named("meanChiSq") = result.meanChiSq,
        Rcpp::Named("diagnostics") = Rcpp::List::create()
    );
}

// these are the functions exposed to the R package

// [[Rcpp::export]]
Rcpp::List cogaps_cpp_from_file(const std::string &data,
const Rcpp::List &allParams,
const std::string &uncertainty,
const Rcpp::NumericVector &indices=Rcpp::NumericVector(1),
const Rcpp::NumericMatrix &fixedMatrix=Rcpp::NumericMatrix(1,1))
{
    return cogapsRun(data, allParams, uncertainty, indices, fixedMatrix);
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &data,
const Rcpp::List &allParams,
const Rcpp::NumericMatrix &uncertainty,
const Rcpp::NumericVector &indices=Rcpp::NumericVector(1),
const Rcpp::NumericMatrix &fixedMatrix=Rcpp::NumericMatrix(1,1))
{
    return cogapsRun(convertRMatrix(data), allParams,
        convertRMatrix(uncertainty), indices, fixedMatrix);
}

// [[Rcpp::export]]
std::string getBuildReport_cpp()
{
    return buildReport();
}
