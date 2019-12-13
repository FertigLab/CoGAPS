#include "GapsParameters.h"
#include "GapsResult.h"
#include "GapsRunner.h"
#include "data_structures/Matrix.h"
#include "file_parser/FileParser.h"
#include "math/Random.h"
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

// convert std::vector of Matrix types to an R list
template <class GenericMatrix>
static Rcpp::List createListOfRMatrices(const std::vector<GenericMatrix> &cppMatrices)
{
    Rcpp::List rMatrices;
    for (unsigned i = 0; i < cppMatrices.size(); ++i)
    {
        rMatrices.push_back(createRMatrix(cppMatrices[i]));
    }
    return rMatrices;
}

////////// converts R parameters to single GapsParameters struct ///////////////

template <class DataType>
GapsParameters getGapsParameters(const DataType &data, const Rcpp::List &allParams)
{
    const Rcpp::S4 &gapsParams(allParams["gaps"]);

    // check if subsetting data
    unsigned subsetDim = Rcpp::as<unsigned>(gapsParams.slot("subsetDim"));
    bool subsetGenes = (subsetDim == 1);
    std::vector<unsigned> subset;
    if (subsetDim > 0)
    {
        Rcpp::IntegerVector subsetR = gapsParams.slot("subsetIndices");
        subset = Rcpp::as< std::vector<unsigned> >(subsetR);
    }

    // create standard CoGAPS parameters struct
    GapsParameters params(data, Rcpp::as<bool>(allParams["transposeData"]),
        subsetDim > 0, subsetGenes, subset);
    params.runningDistributed = subsetDim > 0;
    params.printThreadUsage = !params.runningDistributed;

    // get configuration parameters
    params.maxThreads = Rcpp::as<int>(allParams["nThreads"]);
    params.workerID = Rcpp::as<int>(allParams["workerID"]);
    params.printMessages = Rcpp::as<bool>(allParams["messages"]) && (params.workerID == 1);
    params.outputFrequency = Rcpp::as<int>(allParams["outputFrequency"]);
    params.checkpointOutFile = Rcpp::as<std::string>(allParams["checkpointOutFile"]);
    params.checkpointInterval = Rcpp::as<int>(allParams["checkpointInterval"]);
    params.takePumpSamples = Rcpp::as<bool>(gapsParams.slot("takePumpSamples"));

    // extract model specific parameters from list
    params.seed = Rcpp::as<int>(gapsParams.slot("seed"));
    params.nPatterns = Rcpp::as<int>(gapsParams.slot("nPatterns"));
    params.nIterations = gapsParams.slot("nIterations");
    params.alphaA = Rcpp::as<float>(gapsParams.slot("alphaA"));
    params.alphaP = Rcpp::as<float>(gapsParams.slot("alphaP"));
    params.maxGibbsMassA = Rcpp::as<float>(gapsParams.slot("maxGibbsMassA"));
    params.maxGibbsMassP = Rcpp::as<float>(gapsParams.slot("maxGibbsMassP"));
    params.singleCell = Rcpp::as<bool>(gapsParams.slot("singleCell"));
    params.useSparseOptimization = Rcpp::as<bool>(gapsParams.slot("sparseOptimization"));
    params.asynchronousUpdates = Rcpp::as<bool>(allParams["asynchronousUpdates"]);

    // calculate snapshot frequency
    int nSnapshots = Rcpp::as<int>(allParams["nSnapshots"]);
    if (nSnapshots > 0)
    {
        params.snapshotFrequency = params.nIterations / nSnapshots;
    }

    // determine which phase to take snapshots in
    std::string phase = Rcpp::as<std::string>(allParams["snapshotPhase"]);
    if (phase == "equilibration")
    {
        params.snapshotPhase = GAPS_EQUILIBRATION_PHASE;
    }
    else if (phase == "sampling")
    {
        params.snapshotPhase = GAPS_SAMPLING_PHASE;
    }

    // check if using fixed matrix
    params.whichMatrixFixed = Rcpp::as<char>(gapsParams.slot("whichMatrixFixed"));
    if (params.whichMatrixFixed != 'N')
    {
        params.useFixedPatterns = true;
        Rcpp::NumericMatrix fixedMatrixR = gapsParams.slot("fixedPatterns");
        params.fixedPatterns = convertRMatrix(fixedMatrixR);
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
const DataType &uncertainty)
{
    // convert R parameters to GapsParameters struct
    GapsParameters params(getGapsParameters(data, allParams));
#ifdef GAPS_DEBUG
    params.print();
#endif

    // create GapsRunner, note we must first initialize the random generator
    GapsRandomState randState(params.seed);
    GapsResult result(gaps::run(data, params, uncertainty, &randState));

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
            Rcpp::Named("atomsP") = Rcpp::wrap(result.atomHistoryP),
            Rcpp::Named("pumpStat") = createRMatrix(result.pumpMatrix),
            Rcpp::Named("meanPatternAssignment") = createRMatrix(result.meanPatternAssignment),
            Rcpp::Named("averageQueueLengthA") = result.averageQueueLengthA,
            Rcpp::Named("averageQueueLengthP") = result.averageQueueLengthP,
            Rcpp::Named("totalUpdates") = result.totalUpdates,
            Rcpp::Named("totalRunningTime") = result.totalRunningTime,
            Rcpp::Named("equilibrationSnapshotsA") = createListOfRMatrices(result.equilibrationSnapshotsA),
            Rcpp::Named("equilibrationSnapshotsP") = createListOfRMatrices(result.equilibrationSnapshotsP),
            Rcpp::Named("samplingSnapshotsA") = createListOfRMatrices(result.samplingSnapshotsA),
            Rcpp::Named("samplingSnapshotsP") = createListOfRMatrices(result.samplingSnapshotsP)
        )
    );
}

/////////////////// functions exposed to the R package /////////////////////////

// [[Rcpp::export]]
Rcpp::List cogaps_from_file_cpp(const Rcpp::CharacterVector &data,
const Rcpp::List &allParams,
const Rcpp::Nullable<Rcpp::CharacterVector> &uncertainty)
{
    std::string unc;
    if (uncertainty.isNotNull())
    {
        unc = Rcpp::as<std::string>(Rcpp::CharacterVector(uncertainty));
    }
    return cogapsRun(Rcpp::as<std::string>(data), allParams, unc);
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &data,
const Rcpp::List &allParams,
const Rcpp::Nullable<Rcpp::NumericMatrix> &uncertainty)
{
    Matrix unc;
    if (uncertainty.isNotNull())
    {
        unc = convertRMatrix(Rcpp::NumericMatrix(uncertainty));
    }
    return cogapsRun(convertRMatrix(data), allParams, unc);
}

// [[Rcpp::export]]
std::string getBuildReport_cpp()
{
    return buildReport();
}

// [[Rcpp::export]]
bool checkpointsEnabled_cpp()
{
#ifdef GAPS_DISABLE_CHECKPOINTS
    return false;
#else
    return true;
#endif
}

// [[Rcpp::export]]
bool compiledWithOpenMPSupport_cpp()
{
#ifdef __GAPS_OPENMP__
    return true;
#else
    return false;
#endif
}

// [[Rcpp::export]]
Rcpp::List getFileInfo_cpp(const std::string &path)
{
    FileParser fp(path);
    Rcpp::NumericVector dim(2);
    dim[0] = fp.nRow();
    dim[1] = fp.nCol();
    return Rcpp::List::create(
        Rcpp::Named("dimensions") = dim,
        Rcpp::Named("rowNames") = Rcpp::wrap(fp.rowNames()),
        Rcpp::Named("colNames") = Rcpp::wrap(fp.colNames())
    );
}