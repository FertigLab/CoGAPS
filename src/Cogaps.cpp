#include "GapsRunner.h"
#include "math/SIMD.h"

#include <Rcpp.h>
#include <string>
#include <sstream>

// these are helper functions for converting matrix/vector types
// to and from R objects

static Matrix convertRMatrix(const Rcpp::NumericMatrix &rmat, bool transpose=false)
{
    unsigned nr = transpose ? rmat.ncol() : rmat.nrow();
    unsigned nc = transpose ? rmat.nrow() : rmat.ncol();
    Matrix mat(nr, nc);
    for (unsigned i = 0; i < nr; ++i)
    {
        for (unsigned j = 0; j < nc; ++j)
        {
            mat(i,j) = transpose ? rmat(j,i) : rmat(i,j);
        }
    }
    return mat;
}

template <class Matrix>
static Rcpp::NumericMatrix createRMatrix(const Matrix &mat, bool transpose=false)
{
    unsigned nr = transpose ? mat.nCol() : mat.nRow();
    unsigned nc = transpose ? mat.nRow() : mat.nCol();
    Rcpp::NumericMatrix rmat(nr, nc);
    for (unsigned i = 0; i < nr; ++i)
    {
        for (unsigned j = 0; j < nc; ++j)
        {
            rmat(i,j) = transpose ? mat(j,i) : mat(i,j);
        }
    }
    return rmat;
}

// these helper functions provide an abtracted way for communicating which
// parameters are null between R and C++

static bool isNull(const std::string &file)
{
    return file.empty();
}

static bool isNull(const Matrix &mat)
{
    return mat.nRow() == 1 && mat.nCol() == 1;
}

// needed to create proper size of GapsRunner
unsigned getNumPatterns(const Rcpp::List &allParams)
{
    const Rcpp::S4 &gapsParams(allParams["gaps"]);
    unsigned nPatterns = gapsParams.slot("nPatterns");
    if (!Rf_isNull(allParams["checkpointInFile"]))
    {
        std::string file(Rcpp::as<std::string>(allParams["checkpointInFile"]));
        Archive ar(file, ARCHIVE_READ);
        GapsRng::load(ar);
        ar >> nPatterns;
        ar.close();
    }
    return nPatterns;
}

std::vector<unsigned> getSubsetIndices(const Rcpp::Nullable<Rcpp::IntegerVector> &indices)
{
    if (indices.isNotNull())
    {
        return Rcpp::as< std::vector<unsigned> >(Rcpp::IntegerVector(indices));
    }
    return std::vector<unsigned>(1); // interpreted as null, i.e. will be ignored
}

// return if running distributed, and if so, are we partitioning rows/cols
std::pair<bool, bool> processDistributedParameters(const Rcpp::List &allParams)
{
    const Rcpp::S4 &gapsParams(allParams["gaps"]);
    if (!Rf_isNull(gapsParams.slot("distributed")))
    {
        std::string d(Rcpp::as<std::string>(gapsParams.slot("distributed")));
        GAPS_ASSERT(d == "genome-wide" || d == "single-cell");
        return std::pair<bool, bool>(true, d == "genome-wide");
    }
    return std::pair<bool, bool>(false, false);
}

// this is the main function that creates a GapsRunner and runs CoGAPS

template <class DataType>
static Rcpp::List cogapsRun(const DataType &data, const Rcpp::List &allParams,
const DataType &uncertainty, const Rcpp::Nullable<Rcpp::IntegerVector> &indices,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix, bool isMaster)
{
    // calculate essential parameters needed for constructing GapsRunner
    const Rcpp::S4 &gapsParams(allParams["gaps"]);
    GapsRng::setSeed(gapsParams.slot("seed"));
    unsigned nPatterns = getNumPatterns(allParams); // TODO clarify this sets the checkpoint seed as well
    bool printThreads = !processDistributedParameters(allParams).first;
    bool partitionRows = processDistributedParameters(allParams).second;
    std::vector<unsigned> cIndices(getSubsetIndices(indices));

    // construct GapsRunner
    GapsRunner runner(data, allParams["transposeData"], nPatterns,
        partitionRows, cIndices);

    // set uncertainty
    if (!isNull(uncertainty))
    {
        runner.setUncertainty(uncertainty, allParams["transposeData"],
            partitionRows, cIndices);
    }
    
    // populate GapsRunner from checkpoint file
    if (!Rf_isNull(allParams["checkpointInFile"]))
    {
        std::string file(Rcpp::as<std::string>(allParams["checkpointInFile"]));
        Archive ar(file, ARCHIVE_READ);
        GapsRng::load(ar);
        ar >> runner;
        ar.close();
    }
    else // no checkpoint, populate from given parameters
    {
        // set fixed matrix
        if (fixedMatrix.isNotNull())
        {
            GAPS_ASSERT(!Rf_isNull(allParams["whichMatrixFixed"]));
            std::string which = Rcpp::as<std::string>(allParams["whichMatrixFixed"]);
            runner.setFixedMatrix(which[0], convertRMatrix(Rcpp::NumericMatrix(fixedMatrix), which[0]=='P'));
        }

        // set parameters that would be saved in the checkpoint
        runner.recordSeed(gapsParams.slot("seed"));
        runner.setMaxIterations(gapsParams.slot("nIterations"));
        runner.setSparsity(gapsParams.slot("alphaA"),
            gapsParams.slot("alphaP"), gapsParams.slot("singleCell"));
        runner.setMaxGibbsMass(gapsParams.slot("maxGibbsMassA"),
            gapsParams.slot("maxGibbsMassP"));
    }

    // set parameters that aren't saved in the checkpoint
    runner.setMaxThreads(allParams["nThreads"]);
    runner.setPrintMessages(allParams["messages"] && isMaster);
    runner.setOutputFrequency(allParams["outputFrequency"]);
    runner.setCheckpointOutFile(allParams["checkpointOutFile"]);
    runner.setCheckpointInterval(allParams["checkpointInterval"]);

    // run cogaps
    GapsResult result(runner.run(printThreads));
    
    // write result to file if requested
    if (allParams["outputToFile"] != R_NilValue)
    {
        result.writeToFile(Rcpp::as<std::string>(allParams["outputToFile"]));
    }

    // return R list
    return Rcpp::List::create(
        Rcpp::Named("Amean") = createRMatrix(result.Amean),
        Rcpp::Named("Pmean") = createRMatrix(result.Pmean, true),
        Rcpp::Named("Asd") = createRMatrix(result.Asd),
        Rcpp::Named("Psd") = createRMatrix(result.Psd, true),
        Rcpp::Named("seed") = runner.getSeed(),
        Rcpp::Named("meanChiSq") = result.meanChiSq,
        Rcpp::Named("geneNames") = allParams["geneNames"],
        Rcpp::Named("sampleNames") = allParams["sampleNames"],
        Rcpp::Named("diagnostics") = Rcpp::List::create()
    );
}

// these are the functions exposed to the R package

// [[Rcpp::export]]
Rcpp::List cogaps_cpp_from_file(const Rcpp::CharacterVector &data,
const Rcpp::List &allParams,
const Rcpp::Nullable<Rcpp::CharacterVector> &uncertainty=R_NilValue,
const Rcpp::Nullable<Rcpp::IntegerVector> &indices=R_NilValue,
const Rcpp::Nullable<Rcpp::NumericMatrix> &fixedMatrix=R_NilValue,
bool isMaster=true)
{
    std::string unc = ""; // interpreted as null, i.e. will be ignored
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
    Matrix unc(1,1); // interpreted as null, i.e. will be ignored
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
