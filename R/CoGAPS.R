#' @include class-CogapsParams.R
NULL

#' CoGAPS
#' @name CoGAPS Matrix Factorization Algorithm
#' @docType methods
#' @rdname CoGAPS-methods
#' @description calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix
#' @details Currently, raw count matrices are the only supported R object. For
#'  file types CoGAPS supports csv, tsv, and mtx
#' @param data File name or R object (see details for supported types)
#' @param params CogapsParams object
#' @param uncertainty uncertainty matrix (same supported types as data)
#' @param fixedMatrix data for fixing the values of either the A or P matrix;
#'  used in conjuction with whichMatrixFixed (see CogapsParams)
#' @param checkpointFile name of the checkpoint file
#' @param ... keeps backwards compatibility with arguments from older versions
#' @return CogapsResult object
#' @examples
#' # Running from R object
#' data(GIST)
#' resultA <- CoGAPS(GIST.D)
#' # Running from file name
#' gist_path <- system.file("extdata/GIST.mtx", package="CoGAPS")
#' resultB <- CoGAPS(gist_path)
#' Setting Parameters
#' params <- new("CogapsParams")
#' params <- setParam(params, "nPatterns", 5)
#' resultC <- CoGAPS(GIST.D, params)
#' @importFrom methods new
#' @export
setGeneric("CoGAPS", function(data, params=new("CogapsParams"),
uncertainty=NULL, fixedMatrix=NULL, checkpointFile=NULL, ...)
{
    # parse parameters from ...
    params <- parseOldParams(params, list(...))
    params <- parseDirectParams(params, list(...))

    # call method
    gapsReturnList <- standardGeneric("CoGAPS")

    # convert list to CogapsResult object
    return(CogapsResult(
        Amean       = gapsReturnList$Amean,
        Asd         = gapsReturnList$Asd,
        Pmean       = gapsReturnList$Pmean,
        Psd         = gapsReturnList$Psd,
        seed        = gapsReturnList$seed,
        meanChiSq   = gapsReturnList$meanChiSq,
        diagnostics = gapsReturnList$diagnostics
    )) 
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
#' @importFrom tools file_ext
setMethod("CoGAPS", signature(data="character", params="CogapsParams"),
function(data, params, uncertainty, fixedMatrix, checkpointFile, ...)
{
    # check file extension
    if (!(file_ext(data) %in% c("tsv", "csv", "mtx")))
        stop("unsupported file extension for data")

    # check uncertainty matrix
    if (!is.null(uncertainty))
    {
        if (class(uncertainty) != "character")
            stop("uncertainty must be same data type as data (file name)")
        if (!(file_ext(uncertainty) %in% c("tsv", "csv", "mtx")))
            stop("unsupported file extension for uncertainty")
    }

    # call C++ function
    cogaps_cpp_from_file(data, uncertainty, params@nPatterns,
        params@nIterations, params@outputFrequency, params@seed, params@alphaA,
        params@alphaP, params@maxGibbsMassA, params@maxGibbsMassP,
        params@messages, params@singleCell, params@checkpointOutFile,
        params@nCores)
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod("CoGAPS", signature(data="matrix", params="CogapsParams"),
function(data, params, uncertainty, fixedMatrix, checkpointFile, ...)
{
    # check matrix
    if (!is.null(uncertainty) & class(uncertainty) != "matrix")
        stop("uncertainty must be same data type as data (matrix)")
    checkDataMatrix(data, uncertainty)

    # call C++ function
    cogaps_cpp(data, uncertainty, params@nPatterns,
        params@nIterations, params@outputFrequency, params@seed, params@alphaA,
        params@alphaP, params@maxGibbsMassA, params@maxGibbsMassP,
        params@messages, params@singleCell, params@checkpointOutFile,
        params@nCores)
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod("CoGAPS", signature(data="data.frame", params="CogapsParams"),
function(data, params, uncertainty, fixedMatrix, checkpointFile, ...)
{
    # check matrix
    if (!is.null(uncertainty) & class(uncertainty) != "matrix")
        stop("uncertainty must be matrix when data is a data.frame")
    checkDataMatrix(data.matrix(data), uncertainty)

    # call C++ function
    cogaps_cpp(data.matrix(data), uncertainty, params@nPatterns,
        params@nIterations, params@outputFrequency, params@seed, params@alphaA,
        params@alphaP, params@maxGibbsMassA, params@maxGibbsMassP,
        params@messages, params@singleCell, params@checkpointOutFile,
        params@nCores)
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("CoGAPS", signature(data="SummarizedExperiment", params="CogapsParams"),
function(data, params, uncertainty, fixedMatrix, checkpointFile, ...)
{
    # extract count matrix
    countMatrix = assay(data, "counts")

    # check matrix
    if (!is.null(uncertainty) & class(uncertainty) != "matrix")
        stop("uncertainty must be matrix when data is a SummarizedExperiment")
    checkDataMatrix(countMatrix, uncertainty)

    # call C++ function
    cogaps_cpp(countMatrix, uncertainty, params@nPatterns,
        params@nIterations, params@outputFrequency, params@seed, params@alphaA,
        params@alphaP, params@maxGibbsMassA, params@maxGibbsMassP,
        params@messages, params@singleCell, params@checkpointOutFile,
        params@nCores)
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("CoGAPS", signature(data="SingleCellExperiment", params="CogapsParams"),
function(data, params, uncertainty, fixedMatrix, checkpointFile, ...)
{
    # extract count matrix
    countMatrix = assay(data, "counts")

    # check matrix
    if (!is.null(uncertainty) & class(uncertainty) != "matrix")
        stop("uncertainty must be matrix when data is a SingleCellExperiment")
    checkDataMatrix(countMatrix, uncertainty)

    # call C++ function
    cogaps_cpp(countMatrix, uncertainty, params@nPatterns,
        params@nIterations, params@outputFrequency, params@seed, params@alphaA,
        params@alphaP, params@maxGibbsMassA, params@maxGibbsMassP,
        params@messages, params@singleCell, params@checkpointOutFile,
        params@nCores)
})

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