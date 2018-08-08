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
#' resultA <- CoGAPS(GIST.data_frame)
#'
#' # Running from file name
#' gist_path <- system.file("extdata/GIST.mtx", package="CoGAPS")
#' resultB <- CoGAPS(gist_path)
#'
#' # Setting Parameters
#' params <- new("CogapsParams")
#' params <- setParam(params, "nPatterns", 5)
#' resultC <- CoGAPS(GIST.data_frame, params)
#' @importFrom methods new is
#' @importFrom SummarizedExperiment assay
#' @importFrom utils packageVersion
CoGAPS <- function(data, params=new("CogapsParams"), nThreads=1,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, BPPARAM=NULL,
geneNames=NULL, sampleNames=NULL, matchedPatterns=NULL,
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
        "bpBackend"=BPPARAM,
        "matchedPatterns"=matchedPatterns,
        "outputToFile"=outputToFile,
        "whichMatrixFixed"=NULL # internal parameter
    )
    allParams <- parseExtraParams(allParams, list(...))

    # display start up message for the user
    startupMessage(data, allParams)

    # check file extension
    if (is(data, "character") & !supported(data))
        stop("unsupported file extension for data")

    # check uncertainty matrix
    if (is(data, "character") & !is.null(uncertainty) & !is(uncertainty, "character"))
        stop("uncertainty must be same data type as data (file name)")
    if (is(uncertainty, "character") & !supported(uncertainty))
        stop("unsupported file extension for uncertainty")
    if (!is(data, "character") & !is.null(uncertainty) & !is(uncertainty, "matrix"))
        stop("uncertainty must be a matrix unless data is a file path")
    if (!is(data, "character"))
        checkDataMatrix(data, uncertainty, allParams$gaps)

    # check single cell parameter
    if (!is.null(allParams$gaps@distributed))
        if (allParams$gaps@distributed == "single-cell" & !allParams$gaps@singleCell)
            warning("running single-cell CoGAPS with singleCell=FALSE")

    if (!is.null(allParams$gaps@distributed) & allParams$nThreads > 1)
        stop("can't run multi-threaded and distributed CoGAPS at the same time")

    # convert data to matrix
    if (is(data, "matrix"))
        data <- data
    if (is(data, "data.frame"))
        data <- data.matrix(data)
    else if (is(data, "SummarizedExperiment"))
        data <- SummarizedExperiment::assay(data, "counts")
    else if (is(data, "SingleCellExperiment"))
        data <- SummarizedExperiment::assay(data, "counts")
   
    # determine which function to call cogaps algorithm
    if (!is.null(allParams$gaps@distributed))
        dispatchFunc <- distributedCogaps # genome-wide or single-cell cogaps
    else if (is(data, "character"))
        dispatchFunc <- cogaps_cpp_from_file # data is a file path
    else
        dispatchFunc <- cogaps_cpp # default

    # check if we're running from a checkpoint
    if (!is.null(allParams$checkpointInFile))
    {
        if (!is.null(allParams$gaps@distributed))
            stop("checkpoints not supported for distributed cogaps")
        else
            cat("Running CoGAPS from a checkpoint")
    }

    # get gene/sample names
    if (is.null(geneNames)) geneNames <- getGeneNames(data, allParams$transposeData)
    if (is.null(sampleNames)) sampleNames <- getSampleNames(data, allParams$transposeData)

    nGenes <- ifelse(allParams$transposeData, ncolHelper(data), nrowHelper(data))
    nSamples <- ifelse(allParams$transposeData, nrowHelper(data), ncolHelper(data))

    if (length(geneNames) != nGenes)
        stop("incorrect number of gene names given")
    if (length(sampleNames) != nSamples)
        stop("incorrect number of sample names given")

    allParams$geneNames <- geneNames
    allParams$sampleNames <- sampleNames

    # run cogaps
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