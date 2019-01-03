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
#' @param subsetIndices set of indices to use from the data
#' @param subsetDim which dimension (1=rows, 2=cols) to subset
#' @param BPPARAM BiocParallel backend 
#' @param geneNames vector of names of genes in data
#' @param sampleNames vector of names of samples in data
#' @param fixedPatterns fix either 'A' or 'P' matrix to these values, in the
#' context of distributed CoGAPS (GWCoGAPS/scCoGAPS), the first phase is
#' skipped and fixedPatterns is used for all sets - allowing manual pattern
#' matching, as well as fixed runs of standard CoGAPS
#' @param whichMatrixFixed either 'A' or 'P', indicating which matrix is fixed
#' @param takePumpSamples whether or not to take PUMP samples
#' @param outputToFile name of a file to save the output to, will create 4 files
#' of the form "filename_nPatterns_[Amean, Asd, Pmean, Psd].extension"
#' @param workerID if calling CoGAPS in parallel the worker ID can be specified,
#' only worker 1 prints output and each worker outputs when it finishes, this
#' is not neccesary when using the default parallel methods (i.e. distributed
#' CoGAPS) but only when the user is manually calling CoGAPS in parallel
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
#' @importFrom SummarizedExperiment assay
#' @importFrom utils packageVersion
CoGAPS <- function(data, params=new("CogapsParams"), nThreads=1,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, subsetIndices=NULL, subsetDim=0,
BPPARAM=NULL, geneNames=NULL, sampleNames=NULL, fixedPatterns=NULL,
whichMatrixFixed='N', takePumpSamples=FALSE, outputToFile=NULL, workerID=1, ...)
{
    # store all parameters in a list and parse parameters from ...
    validObject(params)
    allParams <- list("gaps"=params,
        "nThreads"=nThreads,
        "messages"=messages,
        "outputFrequency"=outputFrequency,
        "checkpointOutFile"=checkpointOutFile,
        "checkpointInterval"=checkpointInterval,
        "checkpointInFile"=checkpointInFile,
        "transposeData"=transposeData,
        "subsetIndices"=subsetIndices,
        "subsetDim"=subsetDim,
        "BPPARAM"=BPPARAM,
        "fixedPatterns"=fixedPatterns,
        "whichMatrixFixed"=whichMatrixFixed,
        "takePumpSamples"=takePumpSamples,
        "outputToFile"=outputToFile,
        "workerID"=workerID
    )
    allParams <- parseExtraParams(allParams, list(...))

    # convert data if needed
    if (is(data, "data.frame"))
        data <- data.matrix(data)
    else if (is(data, "SummarizedExperiment"))
        data <- SummarizedExperiment::assay(data, "counts")
    else if (is(data, "SingleCellExperiment"))
        data <- SummarizedExperiment::assay(data, "counts")

    # check that inputs are valid, then read the gene/sample names from the data
    checkInputs(data, uncertainty, allParams)
    allParams <- getNamesFromData(data, allParams, geneNames, sampleNames)
   
    # check if we're running from a checkpoint
    if (!is.null(allParams$checkpointInFile))
    {
        gapsCat(allParams, "Running CoGAPS from a checkpoint\n")
    }

    # determine which function to call cogaps algorithm
    dispatchFunc <- cogaps_cpp # default
    if (!is.null(allParams$gaps@distributed))
        dispatchFunc <- distributedCogaps # genome-wide or single-cell cogaps
    else if (is(data, "character"))
        dispatchFunc <- cogaps_cpp_from_file # data is a file path

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
#' data(GIST)
#' params <- new("CogapsParams")
#' params <- setDistributedParams(params, nSets=2)
#' params <- setParam(params, "nIterations", 100)
#' params <- setParam(params, "nPatterns", 3)
#' result <- scCoGAPS(t(GIST.matrix), params, BPPARAM=BiocParallel::SerialParam())
scCoGAPS <- function(data, params=new("CogapsParams"), nThreads=1,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, subsetIndices=NULL, subsetDim=0,
BPPARAM=NULL, geneNames=NULL, sampleNames=NULL, fixedPatterns=NULL,
whichMatrixFixed='N', takePumpSamples=FALSE, outputToFile=NULL, workerID=1, ...)
{
    params@distributed <- "single-cell"
    params@singleCell <- TRUE
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
        subsetIndices=subsetIndices,
        subsetDim=subsetDim,
        BPPARAM=BPPARAM,
        geneNames=geneNames,
        sampleNames=sampleNames,
        fixedPatterns=fixedPatterns,
        whichMatrixFixed=whichMatrixFixed,
        takePumpSamples=takePumpSamples,
        outputToFile=outputToFile,
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
#' data(GIST)
#' params <- new("CogapsParams")
#' params <- setDistributedParams(params, nSets=2)
#' params <- setParam(params, "nIterations", 100)
#' params <- setParam(params, "nPatterns", 3)
#' result <- GWCoGAPS(GIST.matrix, params, BPPARAM=BiocParallel::SerialParam())
GWCoGAPS <- function(data, params=new("CogapsParams"), nThreads=1,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, subsetIndices=NULL, subsetDim=0,
BPPARAM=NULL, geneNames=NULL, sampleNames=NULL, fixedPatterns=NULL,
whichMatrixFixed='N', takePumpSamples=FALSE, outputToFile=NULL, workerID=1, ...)
{
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
        subsetIndices=subsetIndices,
        subsetDim=subsetDim,
        BPPARAM=BPPARAM,
        geneNames=geneNames,
        sampleNames=sampleNames,
        fixedPatterns=fixedPatterns,
        whichMatrixFixed=whichMatrixFixed,
        takePumpSamples=takePumpSamples,
        outputToFile=outputToFile,
        workerID=workerID,
        ...
    )
}