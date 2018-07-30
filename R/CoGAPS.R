#' @include class-CogapsParams.R
NULL

#' Checks if file is supported
#' @param file path to file
#' @return TRUE if file is supported, FALSE if not
#' @importFrom tools file_ext
supported <- function(file)
{
    if (!is(file, "character"))
        return(FALSE)
    return(tools::file_ext(file) %in% c("tsv", "csv", "mtx"))
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
#' @param ... allows for overwriting parameters in params
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
#' @importFrom methods new is
#' @importFrom SummarizedExperiment assay
CoGAPS <- function(data, params=new("CogapsParams"), nThreads=NULL,
messages=TRUE, outputFrequency=500, uncertainty=NULL,
checkpointOutFile="gaps_checkpoint.out", checkpointInterval=1000,
checkpointInFile=NULL, transposeData=FALSE, ...)
{
    # store all parameters in a list and parse parameters from ...
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
    if (is(data, "character") & !supported(data))
        stop("unsupported file extension for data")

    # check uncertainty matrix
    if (is(data, "character") & !is(uncertainty, "character"))
        stop("uncertainty must be same data type as data (file name)")
    if (is(uncertainty, "character") & !supported(uncertainty))
        stop("unsupported file extension for uncertainty")
    if (!is(data, "character") & !is.null(uncertainty) & !is(uncertainty, "matrix"))
        stop("uncertainty must be a matrix unless data is a file path")
    checkDataMatrix(data, uncertainty, params)

    # convert data to matrix
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
            message("Running CoGAPS from a checkpoint")
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

#' Single Cell CoGAPS
#' @export
#'
#' @description wrapper around single-cell distributed algorithm for CoGAPS
#' @inheritParams CoGAPS
#' @return CogapsResult object
#' @importFrom methods new
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

#' Genome Wide CoGAPS
#' @export
#'
#' @description wrapper around genome-wide distributed algorithm for CoGAPS
#' @inheritParams CoGAPS
#' @return CogapsResult object
#' @importFrom methods new
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

#' parse parameters passed through the ... variable
#'
#' @param allParams list of all parameters
#' @param extraParams list of parameters in ...
#' @return allParams with any valid parameters in extraParams added
#' @note will halt with an error if any parameters in extraParams are invalid
#' @importFrom methods slotNames
parseExtraParams <- function(allParams, extraParams)
{
    # parse direct params
    for (s in slotNames(allParams$gaps))
    {
        if (!is.null(extraParams[[s]]))
        {
            allParams$gaps <- setParam(allParams$gaps, s, extraParams[[s]])
            extraParams[[s]] <- NULL
        }
    }

    # check for unrecognized options
    if (length(extraParams) > 0)
        stop(paste("unrecognized argument:", names(extraParams)[1]))

    return(allParams)
}

#' check that provided data is valid
#'
#' @param data data matrix
#' @param uncertainty uncertainty matrix, can be null
#' @param params CogapsParams object
#' @return throws an error if data has problems
checkDataMatrix <- function(data, uncertainty, params)
{
    if (sum(data < 0) > 0 | sum(uncertainty < 0) > 0)
        stop("negative values in data and/or uncertainty matrix")
    if (nrow(data) <= params@nPatterns | ncol(data) <= params@nPatterns)
        stop("nPatterns must be less than dimensions of data")
    if (sum(uncertainty < 1e-5) > 0)
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
