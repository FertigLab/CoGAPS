#include "GapsRunner.h"
#include "utils/GlobalConfig.h"

#include <Rcpp.h>
#include <string>
#include <sstream>

// this file contains the blueprint for creating a wrapper around the C++
// interface used for running CoGAPS. It exposes some functions to R, has a
// method for converting the R parameters to the standard GapsParameters
// struct, and creates a GapsRunner object. The GapsRunner class manages
// all information from a CoGAPS run and is used to set off the run
// and get return data.

////////////////// functions for converting matrix types ///////////////////////

// convert R to C++ data type
static Matrix convertRMatrix(const Rcpp::NumericMatrix &rmat)
{
    Matrix mat(rmat.nrow(), rmat.ncol());
    for (unsigned i = 0; i < rmat.nrow(); ++i)
    {
        for (unsigned j = 0; j < rmat.ncol(); ++j)
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
const Rcpp::List &allParams, bool isMaster,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix,
const Rcpp::Nullable<Rcpp::IntegerVector> &indices)
{
    // Standard CoGAPS parameters struct
    GapsParameters params(data);

    // get configuration parameters
    params.maxThreads = allParams["nThreads"];
    params.printMessages = allParams["messages"] && isMaster;
    params.transposeData = allParams["transposeData"];
    params.outputFrequency = allParams["outputFrequency"];
    params.checkpointOutFile = Rcpp::as<std::string>(allParams["checkpointOutFile"]);
    params.checkpointInterval = allParams["checkpointInterval"];

    // extract model specific parameters from list
    const Rcpp::S4 &gapsParams(allParams["gaps"]);
    params.seed = gapsParams.slot("seed");
    params.nPatterns = gapsParams.slot("nPatterns");
    params.nIterations = gapsParams.slot("nIterations");
    params.alphaA = gapsParams.slot("alphaA");
    params.alphaP = gapsParams.slot("alphaP");
    params.maxGibbsMassA = gapsParams.slot("maxGibbsMassA");
    params.maxGibbsMassP = gapsParams.slot("maxGibbsMassP");
    params.singleCell = gapsParams.slot("singleCell");

    // check if using fixed matrix
    if (fixedMatrix.isNotNull())
    {
        params.useFixedMatrix = true;
        params.fixedMatrix = convertRMatrix(Rcpp::NumericMatrix(fixedMatrix));
    }

    // check if subsetting data
    if (indices.isNotNull())
    {
        params.subsetData = true;
        params.printThreadUsage = false;

        std::string d(Rcpp::as<std::string>(gapsParams.slot("distributed")));
        params.subsetGenes = (d == "genome-wide");
        params.whichFixedMatrix = (d == "genome-wide") ? 'P' : 'A';

        params.dataIndicesSubset =
            Rcpp::as< std::vector<unsigned> >(Rcpp::IntegerVector(indices));
    }

    // check if using checkpoint file, peek at the saved parameters
    if (!Rf_isNull(allParams["checkpointInFile"]))
    {
        params.checkpointFile = Rcpp::as<std::string>(allParams["checkpointInFile"]);
        params.useCheckPoint = true;
        params.peekCheckpoint(params.checkpointFile);
    }

    return params;
}

////////// main function that creates a GapsRunner and runs CoGAPS /////////////

// note uncertainty matrix gets special treatment since it's the same size as
// the data (potentially large), so we want to avoid copying it into the 
// GapsParameters struct temporarily

template <class DataType>
static Rcpp::List cogapsRun(const DataType &data, const Rcpp::List &allParams,
const DataType &uncertainty, const Rcpp::Nullable<Rcpp::IntegerVector> &indices,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix, bool isMaster)
{
    // convert R parameters to GapsParameters struct
    GapsParameters gapsParams(getGapsParameters(data, allParams, isMaster,
        fixedMatrix, indices));

    // create GapsRunner, note we must first initialize the random generator
    GapsRng::setSeed(gapsParams.seed);
    gaps_printf("Loading Data...");
    GapsRunner runner(data, gapsParams);

    // set uncertainty
    if (!uncertainty.empty())
    {
        runner.setUncertainty(uncertainty, gapsParams);
    }
    gaps_printf("Done!\n");
    
    // run cogaps
    GapsResult result(runner.run());

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
        Rcpp::Named("seed") = gapsParams.seed,
        Rcpp::Named("meanChiSq") = result.meanChiSq,
        Rcpp::Named("geneNames") = allParams["geneNames"],
        Rcpp::Named("sampleNames") = allParams["sampleNames"],
        Rcpp::Named("diagnostics") = Rcpp::List::create()
    );
}

/////////////////// functions exposed to the R package /////////////////////////

// [[Rcpp::export]]
Rcpp::List cogaps_cpp_from_file(const Rcpp::CharacterVector &data,
const Rcpp::List &allParams,
const Rcpp::Nullable<Rcpp::CharacterVector> &uncertainty=R_NilValue,
const Rcpp::Nullable<Rcpp::IntegerVector> &indices=R_NilValue,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix=R_NilValue,
bool isMaster=true)
{
    std::string unc;
    if (uncertainty.isNotNull())
    {
        unc = Rcpp::as<std::string>(Rcpp::CharacterVector(uncertainty));
    }

    return cogapsRun(Rcpp::as<std::string>(data), allParams, unc, indices,
        fixedMatrix, isMaster);
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &data,
const Rcpp::List &allParams,
const Rcpp::Nullable<Rcpp::NumericMatrix> &uncertainty=R_NilValue,
const Rcpp::Nullable<Rcpp::IntegerVector> &indices=R_NilValue,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix=R_NilValue,
bool isMaster=true)
{
    Matrix unc;
    if (uncertainty.isNotNull())
    {
        unc = convertRMatrix(Rcpp::NumericMatrix(uncertainty));
    }

    return cogapsRun(convertRMatrix(data), allParams, unc, indices,
        fixedMatrix, isMaster);
}

// [[Rcpp::export]]
std::string getBuildReport_cpp()
{
    return buildReport();
}
