#include "GapsRunner.h"
#include "utils/GlobalConfig.h"

#include <Rcpp.h>
#include <sstream>
#include <string>

// this file contains the blueprint for creating a wrapper around the C++
// interface used for running CoGAPS. It exposes some functions to R, has a
// method for converting the R parameters to the standard GapsParameters
// and calls gaps::run

////////////////// functions for converting matrix types ///////////////////////

// convert R to C++ data type
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

// convert C++ to R data type
template <class GenericMatrix>
static Rcpp::NumericMatrix createRMatrix(const GenericMatrix &mat)
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

////////// converts R parameters to single GapsParameters struct ///////////////

template <class DataType>
GapsParameters getGapsParameters(const DataType &data,
const Rcpp::List &allParams, unsigned workerID,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix,
const Rcpp::Nullable<Rcpp::IntegerVector> &indices)
{
    // check if subsetting data
    const Rcpp::S4 &gapsParams(allParams["gaps"]);
    bool subsetGenes = false;
    std::vector<unsigned> subset;
    if (indices.isNotNull())
    {
        std::string d(Rcpp::as<std::string>(gapsParams.slot("distributed")));
        subsetGenes = (d == "genome-wide");
        subset = Rcpp::as< std::vector<unsigned> >(Rcpp::IntegerVector(indices));
    }

    // create standard CoGAPS parameters struct
    GapsParameters params(data, allParams["transposeData"], indices.isNotNull(),
        subsetGenes, subset);
    params.printThreadUsage = !indices.isNotNull();
    params.whichFixedMatrix = indices.isNotNull() ? (subsetGenes ? 'P' : 'A') : 'N';

    // get configuration parameters
    params.maxThreads = allParams["nThreads"];
    params.printMessages = allParams["messages"] && (workerID == 1);
    params.workerID = workerID;
    params.outputFrequency = allParams["outputFrequency"];
    params.checkpointOutFile = Rcpp::as<std::string>(allParams["checkpointOutFile"]);
    params.checkpointInterval = allParams["checkpointInterval"];

    // extract model specific parameters from list
    params.seed = gapsParams.slot("seed");
    params.nPatterns = gapsParams.slot("nPatterns");
    params.nIterations = gapsParams.slot("nIterations");
    params.alphaA = gapsParams.slot("alphaA");
    params.alphaP = gapsParams.slot("alphaP");
    params.maxGibbsMassA = gapsParams.slot("maxGibbsMassA");
    params.maxGibbsMassP = gapsParams.slot("maxGibbsMassP");
    params.singleCell = gapsParams.slot("singleCell");
    params.useSparseOptimization = gapsParams.slot("sparseOptimization");

    // check if using fixed matrix
    if (fixedMatrix.isNotNull())
    {
        params.useFixedMatrix = true;
        params.fixedMatrix = convertRMatrix(Rcpp::NumericMatrix(fixedMatrix));
    }

    // check if using checkpoint file, peek at the saved parameters
    if (!Rf_isNull(allParams["checkpointInFile"]))
    {
        params.checkpointFile = Rcpp::as<std::string>(allParams["checkpointInFile"]);
        params.useCheckPoint = true;
    }

    return params;
}

////////////////////// main function that runs CoGAPS //////////////////////////

// note uncertainty matrix gets special treatment since it's the same size as
// the data (potentially large), so we want to avoid copying it into the 
// GapsParameters struct temporarily

template <class DataType>
static Rcpp::List cogapsRun(const DataType &data, const Rcpp::List &allParams,
const DataType &uncertainty, const Rcpp::Nullable<Rcpp::IntegerVector> &indices,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix, unsigned workerID)
{
    // convert R parameters to GapsParameters struct
    GapsParameters params(getGapsParameters(data, allParams, workerID,
        fixedMatrix, indices));

    // create GapsRunner, note we must first initialize the random generator
    GapsRandomState randState(params.seed);
    GapsResult result(gaps::run(data, params, uncertainty, &randState));

    // write result to file if requested
    if (allParams["outputToFile"] != R_NilValue)
    {
        result.writeToFile(Rcpp::as<std::string>(allParams["outputToFile"]));
    }
    
    // return R list
    return Rcpp::List::create(
        Rcpp::Named("Amean") = createRMatrix(result.Amean),
        Rcpp::Named("Pmean") = createRMatrix(result.Pmean),
        Rcpp::Named("Asd") = createRMatrix(result.Asd),
        Rcpp::Named("Psd") = createRMatrix(result.Psd),
        Rcpp::Named("seed") = params.seed,
        Rcpp::Named("meanChiSq") = result.meanChiSq,
        Rcpp::Named("geneNames") = allParams["geneNames"],
        Rcpp::Named("sampleNames") = allParams["sampleNames"],
        Rcpp::Named("diagnostics") = Rcpp::List::create(
            Rcpp::Named("chisq") = Rcpp::wrap(result.chisqHistory),
            Rcpp::Named("atomsA") = Rcpp::wrap(result.atomHistoryA),
            Rcpp::Named("atomsP") = Rcpp::wrap(result.atomHistoryP)
        )
    );
}

/////////////////// functions exposed to the R package /////////////////////////

// [[Rcpp::export]]
Rcpp::List cogaps_cpp_from_file(const Rcpp::CharacterVector &data,
const Rcpp::List &allParams,
const Rcpp::Nullable<Rcpp::CharacterVector> &uncertainty=R_NilValue,
const Rcpp::Nullable<Rcpp::IntegerVector> &indices=R_NilValue,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix=R_NilValue,
unsigned workerID=1)
{
    std::string unc;
    if (uncertainty.isNotNull())
    {
        unc = Rcpp::as<std::string>(Rcpp::CharacterVector(uncertainty));
    }

    return cogapsRun(Rcpp::as<std::string>(data), allParams, unc, indices,
        fixedMatrix, workerID);
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &data,
const Rcpp::List &allParams,
const Rcpp::Nullable<Rcpp::NumericMatrix> &uncertainty=R_NilValue,
const Rcpp::Nullable<Rcpp::IntegerVector> &indices=R_NilValue,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix=R_NilValue,
unsigned workerID=1)
{
    Matrix unc;
    if (uncertainty.isNotNull())
    {
        unc = convertRMatrix(Rcpp::NumericMatrix(uncertainty));
    }

    return cogapsRun(convertRMatrix(data), allParams, unc, indices,
        fixedMatrix, workerID);
}

// [[Rcpp::export]]
std::string getBuildReport_cpp()
{
    return buildReport();
}
