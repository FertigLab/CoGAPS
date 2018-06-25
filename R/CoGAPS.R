#' @include class-CoGapsParams.R
NULL

#' CoGAPS Matrix Factorization Algorithm
#' @export
#' @docType methods
#' @rdname CoGAPS-methods
#' 
#' @description calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix
#' @details Currently, raw count matrices are the only supported R object. For
#'  file types CoGAPS supports csv, tsv, and mtx
#' @param data File name or R object (see details for supported types)
#' @param params CoGapsParams object
#' @param uncertainty uncertainty matrix (same supported types as data)
#' @param fixedMatrix data for fixing the values of either the A or P matrix;
#'  used in conjuction with whichMatrixFixed (see CoGapsParams)
#' @param checkpointFile name of the checkpoint file
#' @param ... keeps backwards compatibility with arguments from older versions
#' @return CoGapsResult object
#' @importFrom methods new
#' @examples
#' data(SimpSim)
#' params <- new("CoGapsParams")
#' result <- CoGAPS(SimpSim.D, params)
setGeneric("CoGAPS", function(data, params=new("CoGapsParams"),
uncertainty=NULL, fixedMatrix=NULL, checkpointFile=NULL, ...)
{
    # parse parameters from ...
    params <- getTempParams(params, list(...))

    # call method
    standardGeneric("CoGAPS")
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod("CoGAPS", signature(data="matrix", params="CoGapsParams"),
function(data, params, uncertainty, fixedMatrix,
checkpointFile, ...)
{
    # check matrix
    if (class(uncertainty) != "matrix")
        stop("Uncertainty must be same data type as data (matrix)")
    checkDataMatrix(data, uncertainty)

    # call C++ function
    cogaps_cpp(data, uncertainty, nFactor, nIter, nOutputs, seed, alphaA,
        alphaP, maxGibbsMassA, maxGibbsMassP, messages, singleCell,
        checkpointFile))
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
#' @importFrom tools file_ext
setMethod("CoGAPS", signature(data="character", params="CoGapsParams"),
function(data, params, uncertainty, fixedMatrix,
checkpointFile, ...)
{
    # check file extension
    if (!(file_ext(data) %in% c("tsv", "csv", "mtx")))
        stop("unsupported file extension for data")

    # check uncertainty matrix
    if (!is.null(uncertainty))
    {
        if (!(file_ext(uncertainty) %in% c("tsv", "csv", "mtx")))
            stop("unsupported file extension for uncertainty")
        if (class(uncertainty) != "character")
            stop("Uncertainty must be same data type as data (file name)")
    }

    # call C++ function
    cogaps_cpp_from_file(data, uncertainty, nFactor, nIter, nOutputs, seed,
        alphaA, alphaP, maxGibbsMassA, maxGibbsMassP, messages, singleCell,
        checkpointFile))
})

#' Information About Package Compilation
#'
#' @details returns information about how the package was compiled, i.e. which
#'  compiler/version was used, which compile time options were enabled, etc...
#' @return string containing build report
#' @examples
#'  CoGAPS::buildReport()
#' @export
buildReport <- function()
{
    getBuildReport_cpp()
}

#' Display Information About Package Compilation
#'
#' @details displays information about how the package was compiled, i.e. which
#'  compiler/version was used, which compile time options were enabled, etc...
#' @return display builds information
#' @examples
#'  CoGAPS::displayBuildReport()
#' @export
displayBuildReport <- function()
{
    cat(getBuildReport_cpp())
}
