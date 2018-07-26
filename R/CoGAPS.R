#' @include class-CogapsParams.R
NULL

#' CoGAPS Matrix Factorization Algorithm
#' @export 
#' @docType methods
#' @rdname CoGAPS-methods
#'
#' @description calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix
#' @details For file types CoGAPS supports csv, tsv, and mtx
#' @param data File name or R object (see details for supported types)
#' @param params CogapsParams object
#' @param uncertainty uncertainty matrix (same supported types as data)
#' @param fixedMatrix data for fixing the values of either the A or P matrix;
#'  used in conjuction with whichMatrixFixed (see CogapsParams)
#' @param checkpointInFile name of the checkpoint file
#' @param ... keeps backwards compatibility with arguments from older versions
#' @return CogapsResult object
#' @examples
#' # Running from R object
#' data(GIST)
#' resultA <- CoGAPS(GIST.D)
#'
#' # Running from file name
#' gist_path <- system.file("extdata/GIST.mtx", package="CoGAPS")
#' resultB <- CoGAPS(gist_path)
#'
#' Setting Parameters
#' params <- new("CogapsParams")
#' params <- setParam(params, "nPatterns", 5)
#' resultC <- CoGAPS(GIST.D, params)
#' @importFrom methods new
CoGAPS <- function(data, params=new("CogapsParams"), nThreads=NULL,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, ...)
{
    # parse parameters from ...
    allParams <- list("gaps"=params,
        "nThreads"=nThreads,
        "messages"=messages,
        "outputFrequency"=outputFrequency,
        "checkpointOutFile"=checkpointOutFile,
        "checkpointInterval"=checkpointInterval,
        "checkpointInFile"=checkpointInFile,
        "transposeData"=transposeData,
        "whichMatrixFixed"="" # internal parameter
    )
    allParams <- parseExtraParams(allParams, list(...))

    # check file extension
    if (class(data) == "character" & !(file_ext(data) %in% c("tsv", "csv", "mtx")))
        stop("unsupported file extension for data")

    # check uncertainty matrix
    if (class(data) == "character" & class(uncertainty) != "character")
        stop("uncertainty must be same data type as data (file name)")
    if (nchar(uncertainty) > 0 & !(file_ext(uncertainty) %in% c("tsv", "csv", "mtx")))
        stop("unsupported file extension for uncertainty")

    # check matrix
    if (class(uncertainty) != "matrix")
        stop("uncertainty must be same data type as data (matrix)")
    checkDataMatrix(data, uncertainty, params)

    # check if uncertainty is null
    if (is.null(uncertainty) & class(data) == "character")
        uncertainty <- ""
    else if (is.null(uncertainty))
        uncertainty <- matrix(0)

    # convert data to matrix
    if (class(data) == "data.frame")
        data <- data.matrix(data)
    else if (class(data) == "SummarizedExperiment")
        data <- assay(data, "counts")
    else if (class(data) == "SingleCellExperiment")
        data <- assay(data, "counts")

    # determine which function to call cogaps algorithm
    if (!is.null(callParams$gapsParams@distributed))
        dispatchFunc <- distributedCogaps # genome-wide or single-cell cogaps
    else if (class(data) == "character")
        dispatchFunc <- cogaps_cpp_from_file # data is a file name
    else
        dispatchFunc <- cogaps_cpp # default

    # check if we're running from a checkpoint
    if (!is.null(allParams$checkpointInFile))
    {
        if (!is.null(callParams$gapsParams@distributed))
            stop("checkpoints not supported for distributed cogaps")
        else
            message("Running CoGAPS from a checkpoint")
    }
    else
    {
        allParams$checkpointInFile <- ""
    }

    # run cogaps
    gapsReturnList <- dispatchFunc(data, allParams, uncertainty)

    # convert list to CogapsResult object
    return(new("CogapsResult",
        Amean       = gapsReturnList$Amean,
        Asd         = gapsReturnList$Asd,
        Pmean       = gapsReturnList$Pmean,
        Psd         = gapsReturnList$Psd,
        seed        = gapsReturnList$seed,
        meanChiSq   = gapsReturnList$meanChiSq,
        diagnostics = list("diag"=gapsReturnList$diagnostics, "params"=params)
    ))
}

scCoGAPS <- function(data, params=new("CogapsParams"), nThreads=NULL,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, ...)
{
    params@distributed <- "single-cell"
    CoGAPS(data, params, nThreads, messages, outputFrequency, uncertainty,
        checkpointOutFile, checkpointInterval, checkpointInFile, transposeData,
        ...)
}

GWCoGAPS <- function(data, params=new("CogapsParams"), nThreads=NULL,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, ...)
{
    params@distributed <- "genome-wide"
    CoGAPS(data, params, nThreads, messages, outputFrequency, uncertainty,
        checkpointOutFile, checkpointInterval, checkpointInFile, transposeData,
        ...)
}   

parseExtraParams <- function(allParams, extraParams)
{


}

#' Check that provided data is valid
#'
#' @param data data matrix
#' @param uncertainty uncertainty matrix
#' @return throws an error if data has problems
checkDataMatrix <- function(data, uncertainty, params)
{
    if (sum(data < 0) > 0 | sum(uncertainty < 0) > 0)
        stop("negative values in data and/or uncertainty matrix")
    if (nrow(data) <= params@nPatterns | ncol(data) <= params@nPatterns)
        stop("nPatterns must be less than dimensions of data")
    if (sum(dim(uncertainty)) != 2 & sum(uncertainty < 1e-5) > 0)
        warning("small values in uncertainty matrix detected")
}

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
