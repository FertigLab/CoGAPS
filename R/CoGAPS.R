#' @include class-CogapsParams.R
NULL

#' Information About Package Compilation
#' @export
#'
#' @details returns information about how the package was compiled, i.e. which
#'  compiler/version was used, which compile time options were enabled, etc...
#' @return string containing build report
#' @examples
#' CoGAPS::buildReport()
buildReport <- function()
{
    getBuildReport_cpp()
}

#' Check if package was built with checkpoints enabled
#' @export
#'
#' @return true/false if checkpoints are enabled
#' @examples
#' CoGAPS::checkpointsEnabled()
checkpointsEnabled <- function()
{
    checkpointsEnabled_cpp()
}

#' Check if compiler supported OpenMP
#' @export
#'
#' @return true/false if OpenMP was supported
#' @examples
#' CoGAPS::compiledWithOpenMPSupport()
compiledWithOpenMPSupport <- function()
{
    compiledWithOpenMPSupport_cpp()
}

#' CoGAPS Matrix Factorization Algorithm
#' @export 
#'
#' @description calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix
#' @details The supported R types are: matrix, data.frame, SummarizedExperiment,
#' SingleCellExperiment. The supported file types are csv, tsv, and mtx.
#' @param data File name or R object (see details for supported types)
#' @param params CogapsParams object
#' @param nThreads maximum number of threads to run on
#' @param messages T/F for displaying output
#' @param outputFrequency number of iterations between each output (set to 0 to
#' disable status updates, other output is controlled by @code messages)
#' @param uncertainty uncertainty matrix - either a matrix or a supported
#' file type
#' @param checkpointOutFile name of the checkpoint file to create
#' @param checkpointInterval number of iterations between each checkpoint (set
#' to 0 to disable checkpoints)
#' @param checkpointInFile if this is provided, CoGAPS runs from the checkpoint
#' contained in this file
#' @param transposeData T/F for transposing data while reading it in - useful
#' for data that is stored as samples x genes since CoGAPS requires data to be
#' genes x samples
#' @param BPPARAM BiocParallel backend 
#' @param workerID if calling CoGAPS in parallel the worker ID can be specified,
#' only worker 1 prints output and each worker outputs when it finishes, this
#' is not neccesary when using the default parallel methods (i.e. distributed
#' CoGAPS) but only when the user is manually calling CoGAPS in parallel
#' @param asynchronousUpdates enable asynchronous updating which allows for multi-threaded runs
#' @param nSnapshots how many snapshots to take in each phase, setting this to 0 disables
#' snapshots
#' @param snapshotPhase which phase to take snapsjots in e.g. "equilibration", "sampling",
#' "all"
#' @param ... allows for overwriting parameters in params
#' @return CogapsResult object
#' @examples
#' # Running from R object
#' data(GIST)
#' resultA <- CoGAPS(GIST.data_frame, nIterations=25)
#'
#' # Running from file name
#' gist_path <- system.file("extdata/GIST.mtx", package="CoGAPS")
#' resultB <- CoGAPS(gist_path, nIterations=25)
#'
#' # Setting Parameters
#' params <- new("CogapsParams")
#' params <- setParam(params, "nPatterns", 3)
#' resultC <- CoGAPS(GIST.data_frame, params, nIterations=25)
#' @importFrom methods new is
CoGAPS <- function(data, params=new("CogapsParams"), nThreads=1, messages=TRUE,
outputFrequency=1000, uncertainty=NULL, checkpointOutFile="gaps_checkpoint.out",
checkpointInterval=0, checkpointInFile=NULL, transposeData=FALSE,
BPPARAM=NULL, workerID=1, asynchronousUpdates=TRUE, nSnapshots=0,
snapshotPhase='sampling', ...)
{
    # pre-process inputs
    if (is(data, "character"))
        dataName <- data
    else
        dataName <- deparse(substitute(data))
    data <- getValueOrRds(data)
    data <- convertDataToMatrix(data)
    params <- getValueOrRds(params)
    validObject(params)

    # check OpenMP support
    if (!compiledWithOpenMPSupport())
    {
        if (asynchronousUpdates & nThreads > 1)
            warning("requesting multi-threaded version of CoGAPS but compiler did not support OpenMP")
        asynchronousUpdates = FALSE
        nThreads = 1
    }

    # store all parameters in a list and parse parameters from ...
    allParams <- list("gaps"=params,
        "nThreads"=nThreads,
        "messages"=messages,
        "outputFrequency"=outputFrequency,
        "nSnapshots"=nSnapshots,
        "snapshotPhase"=snapshotPhase,
        "checkpointOutFile"=checkpointOutFile,
        "checkpointInterval"=checkpointInterval,
        "checkpointInFile"=checkpointInFile,
        "geneNames"=NULL, # the gene/sample names in the params object are kept
        "sampleNames"=NULL, # as a reference, these are the values actually used
        "transposeData"=transposeData,
        "BPPARAM"=BPPARAM,
        "outputToFile"=NULL,
        "workerID"=workerID,
        "asynchronousUpdates"=asynchronousUpdates,
        "dataName"=dataName
    )
    allParams <- parseExtraParams(allParams, list(...))
    allParams <- getDimNames(data, allParams)
    checkInputs(data, uncertainty, allParams)

    # check if we're running from a checkpoint
    if (!is.null(allParams$checkpointInFile))
    {
        gapsCat(allParams, "Running CoGAPS from a checkpoint\n")
    }

    # determine function to call cogaps algorithm
    dispatchFunc <- cogaps_cpp # default
    if (!is.null(allParams$gaps@distributed))
        dispatchFunc <- distributedCogaps # genome-wide or single-cell cogaps
    else if (is(data, "character"))
        dispatchFunc <- cogaps_from_file_cpp # data is a file path

    # run cogaps
    startupMessage(data, allParams)
    gapsReturnList <- dispatchFunc(data, allParams, uncertainty)
    return(createCogapsResult(gapsReturnList, allParams))
}

#' Single Cell CoGAPS
#' @export
#'
#' @description wrapper around single-cell distributed algorithm for CoGAPS
#' @inheritParams CoGAPS
#' @return CogapsResult object
#' @importFrom methods new
#' @examples
#' \dontrun{
#' data(GIST)
#' params <- new("CogapsParams")
#' params <- setDistributedParams(params, nSets=2)
#' params <- setParam(params, "nIterations", 100)
#' params <- setParam(params, "nPatterns", 3)
#' result <- scCoGAPS(t(GIST.matrix), params, BPPARAM=BiocParallel::SerialParam())
#' }
scCoGAPS <- function(data, params=new("CogapsParams"), nThreads=1, messages=TRUE,
outputFrequency=500, uncertainty=NULL, checkpointOutFile="gaps_checkpoint.out",
checkpointInterval=1000, checkpointInFile=NULL, transposeData=FALSE,
BPPARAM=NULL, workerID=1, asynchronousUpdates=FALSE, ...)
{
    warning(paste("scCoGAPS is deprecated, use the main function CoGAPS",
        "with the argument: distributed=\"single-cell\""))
    params@distributed <- "single-cell"
    CoGAPS(
        data=data,
        params=params,
        nThreads=nThreads,
        messages=messages,
        outputFrequency=outputFrequency,
        uncertainty=uncertainty,
        checkpointOutFile=checkpointOutFile,
        checkpointInterval=checkpointInterval,
        checkpointInFile=checkpointInFile,
        transposeData=transposeData,
        BPPARAM=BPPARAM,
        workerID=workerID,
        ...
    )
}

#' Genome Wide CoGAPS
#' @export
#'
#' @description wrapper around genome-wide distributed algorithm for CoGAPS
#' @inheritParams CoGAPS
#' @return CogapsResult object
#' @importFrom methods new
#' @examples
#' \dontrun{
#' data(GIST)
#' params <- new("CogapsParams")
#' params <- setDistributedParams(params, nSets=2)
#' params <- setParam(params, "nIterations", 100)
#' params <- setParam(params, "nPatterns", 3)
#' result <- GWCoGAPS(GIST.matrix, params, BPPARAM=BiocParallel::SerialParam())
#' }
GWCoGAPS <- function(data, params=new("CogapsParams"), nThreads=1, messages=TRUE,
outputFrequency=500, uncertainty=NULL, checkpointOutFile="gaps_checkpoint.out",
checkpointInterval=1000, checkpointInFile=NULL, transposeData=FALSE,
BPPARAM=NULL, workerID=1, asynchronousUpdates=FALSE, ...)
{
    warning(paste("GWCoGAPS is deprecated, use the main function CoGAPS",
        "with the argument: distributed=\"genome-wide\""))
    params@distributed <- "genome-wide"
    CoGAPS(
        data=data,
        params=params,
        nThreads=nThreads,
        messages=messages,
        outputFrequency=outputFrequency,
        uncertainty=uncertainty,
        checkpointOutFile=checkpointOutFile,
        checkpointInterval=checkpointInterval,
        checkpointInFile=checkpointInFile,
        transposeData=transposeData,
        BPPARAM=BPPARAM,
        workerID=workerID,
        ...
    )
}
