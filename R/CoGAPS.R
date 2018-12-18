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
#' @param BPPARAM BiocParallel backend 
#' @param geneNames vector of names of genes in data
#' @param sampleNames vector of names of samples in data
#' @param matchedPatterns manually matched patterns for distributed CoGAPS
#' @param outputToFile name of a file to save the output to, will create 4 files
#' of the form "filename_nPatterns_[Amean, Asd, Pmean, Psd].extension"
#' @param ... allows for overwriting parameters in params
#' @return CogapsResult object
#' @examples
#' # Running from R object
#' data(GIST)
#' resultA <- CoGAPS(GIST.data_frame, nIterations=100)
#'
#' # Running from file name
#' gist_path <- system.file("extdata/GIST.mtx", package="CoGAPS")
#' resultB <- CoGAPS(gist_path, nIterations=100)
#'
#' # Setting Parameters
#' params <- new("CogapsParams")
#' params <- setParam(params, "nPatterns", 5)
#' resultC <- CoGAPS(GIST.data_frame, params, nIterations=100)
#' @importFrom methods new is
#' @importFrom SummarizedExperiment assay
#' @importFrom utils packageVersion
CoGAPS <- function(data, params=new("CogapsParams"), nThreads=1,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, subsetData=NULL, BPPARAM=NULL,
geneNames=NULL, sampleNames=NULL, fixedPatterns=NULL, whichMatrixFixed='N',
outputToFile=NULL, ...)
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
        "subsetData"=subsetData,
        "bpBackend"=BPPARAM,
        "fixedPatterns"=fixedPatterns,
        "whichMatrixFixed"=whichMatrixFixed,
        "outputToFile"=outputToFile,
    )
    allParams <- parseExtraParams(allParams, list(...))

    # check that inputs are valid, then read the gene/sample names from the data
    checkInputs(data, uncertainty, allParams)
    allParams <- getNamesFromData(data, allParams, geneNames, sampleNames)
   
    # check if we're running from a checkpoint
    if (!is.null(allParams$checkpointInFile))
    {
        if (!is.null(allParams$gaps@distributed))
            stop("checkpoints not supported for distributed cogaps")
        else
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

    # convert list to CogapsResult object
    return(new("CogapsResult",
        Amean       = gapsReturnList$Amean,
        Asd         = gapsReturnList$Asd,
        Pmean       = gapsReturnList$Pmean,
        Psd         = gapsReturnList$Psd,
        meanChiSq   = gapsReturnList$meanChiSq,
        geneNames   = gapsReturnList$geneNames,
        sampleNames = gapsReturnList$sampleNames,
        diagnostics = append(gapsReturnList$diagnostics,
            list("params"=allParams$gaps, "version"=utils::packageVersion("CoGAPS")))
    ))
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
#' result <- scCoGAPS(t(GIST.matrix), BPPARAM=BiocParallel::SerialParam(), nIterations=100)
scCoGAPS <- function(data, params=new("CogapsParams"), nThreads=1,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, BPPARAM=NULL,
geneNames=NULL, sampleNames=NULL, matchedPatterns=NULL, ...)
{
    params@distributed <- "single-cell"
    params@singleCell <- TRUE
    CoGAPS(data, params, nThreads, messages, outputFrequency, uncertainty,
        checkpointOutFile, checkpointInterval, checkpointInFile, transposeData,
        BPPARAM, geneNames, sampleNames, matchedPatterns, ...)
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
#' result <- GWCoGAPS(GIST.matrix, BPPARAM=BiocParallel::SerialParam(), nIterations=100)
GWCoGAPS <- function(data, params=new("CogapsParams"), nThreads=1,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, BPPARAM=NULL,
geneNames=NULL, sampleNames=NULL, matchedPatterns=NULL, ...)
{
    params@distributed <- "genome-wide"
    CoGAPS(data, params, nThreads, messages, outputFrequency, uncertainty,
        checkpointOutFile, checkpointInterval, checkpointInFile, transposeData,
        BPPARAM, geneNames, sampleNames, matchedPatterns, ...)
}