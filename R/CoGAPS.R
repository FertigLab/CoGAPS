#' CoGAPS Matrix Factorization Algorithm
#' @export
#' @docType methods
#' @rdname CoGAPS-methods
#' 
#' @details calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix
#' @param D data
#' @param S uncertainty
#' @param params object of type CoGapsParams
#' @param ... keeps backwards compatibility with arguments from older versions
#' @return object of type CoGapsResult
#' @importFrom methods new
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D, nFactor=3, nOutputs=250)
setGeneric('CoGAPS', function(D, S=NULL, params, ...)
    {standardGeneric('RVsharing')})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(D='matrix', params='CoGapsParams'),
function(D, S=NULL, params, ...)
{
    # default to 10% of signal if no uncertainty matrix passed in by user
    if (is.null(S))
    {
        S <- pmax(0.1 * D, 0.1)
    }

    return(RunCoGAPS(D, S, params))
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(D='data.frame', params='CoGapsParams'),
function(D, S=NULL, params, ...)
{
    # TODO(Hyejune) convert data.frame to matrix and set default value for
    # S matrix if it's null - then call RunCoGAPS
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(D='SummarizedExperiemnt', params='CoGapsParams'),
function(D, S=NULL, params, ...)
{
    # TODO(Hyejune) extract count matrix from SE object and set default value
    # S matrix if it's null - then call RunCoGAPS
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(D='SingleCellExperiemnt', params='CoGapsParams'),
function(D, S=NULL, params, ...)
{
    # TODO(Hyejune) extract count matrix from SCE object and set default value
    # S matrix if it's null - then call RunCoGAPS
})

RunCoGAPS <- function(D, S, params)
{
    thresholdEnum <- c("unique", "cut")
    result <- cogaps_cpp(D, S, params@nFactor, params@nIter, params@nIter / 10,
        params@nIter, params@outputFrequency, params@alphaA, params@alphaP,
        params@maxGibbmassA, params@maxGibbmassP, params@seed, params@messages,
        params@singleCellRNASeq, params@whichMatrixFixed, params@fixedMatrix,
        params@checkpointInterval, params@checkpointFile, 0, 0, params@nCores)

    return(new('CoGapsResult', result$Amean, result$Asd, result$Pmean,
        result$Psd))
}

# TODO(Tom) handle instance where function is called with deprecated parametres