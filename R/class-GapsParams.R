#' @title GapsParams
#' @description Parameters for running CoGAPS
#'
#' @slot nFactor number of patterns (basis vectors, metagenes)
#' @slot nEquil number of iterations for burn-in
#' @slot nSample number of iterations for sampling
#' @slot nOutput how often (number of iterations) to print status in R console
#' @slot numSnapshots the number of individual samples to capture
#' @slot alphaA sparsity parameter for A domain
#' @slot alphaP sparsity parameter for P domain
#' @slot maxGibbmassA limit truncated normal to max size
#' @slot maxGibbmassP limit truncated normal to max size
#' @slot seed positive values are kept, negative values will be overwritten
#'  by a seed generated from the current time
#' @slot messages display messages during the run
#' @slot singleCelLRNASeq indicate if the data is single cell RNA-seq data
#' @slot fixedPatterns matrix of fixed values in either A or P matrix
#' @slot whichMatrixFixed characted indicating whether A or P contains
#'  the fixed patterns
#' @slot checkpointInterval how often (number of seconds) to create a checkpoint file
#' @slot nCores max number of cores to run on
#' @export
setClass('GapsParams', slots = c(
    nFactor = 'integer',
    nEquil = 'integer',
    nEquilCool = 'integer',
    nSample = 'integer',
    nOutput = 'integer',
    nSnapshots = 'integer',
    alphaA = 'numeric',
    alphaP = 'numeric',
    maxGibbmassA = 'numeric',
    maxGibbmassP = 'numeric',
    seed = 'integer',
    messages  = 'logical',
    singleCellRNASeq = 'logical',
    fixedPatterns = 'matrix',
    whichMatrixFixed = 'character',
    checkpointInterval = 'integer',
    nCores = 'integer'
))

#' @importFrom methods callNextMethod
setMethod('initialize', 'GapsParams', 
    function(.Object, nFactor, nEquil, nSample, ...)
    {
        .Object@nFactor <- as.integer(nFactor)
        .Object@nEquil <- as.integer(nEquil)
        .Object@nEquil <- as.integer(floor(nEquil/10))
        .Object@nSample <- as.integer(nSample)
        .Object@nOutput <- as.integer(1000)
        .Object@nSnapshots <- as.integer(0)
        .Object@alphaA <- 0.01
        .Object@alphaP <- 0.01
        .Object@maxGibbmassA <- 100.0
        .Object@maxGibbmassP <- 100.0
        .Object@seed <- as.integer(-1)
        .Object@messages <- TRUE
        .Object@singleCellRNASeq <- FALSE
        .Object@fixedPatterns <- matrix(0)
        .Object@whichMatrixFixed <- 'N'
        .Object@checkpointInterval <- as.integer(0)
        .Object@nCores <- as.integer(1)

        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

setValidity('GapsParams',
    function(object)
    {
        if (object@whichMatrixFixed == 'P' & nrow(object@fixedPatterns) > object@nFactor)
            'number of fixed patterns greater than nFactor'
        if (object@whichMatrixFixed == 'A' & ncol(object@fixedPatterns) > object@nFactor)
            'number of fixed patterns greater than nFactor'
    }
)

##################### Generics ###################

#' @export
setGeneric('getParam', function(object, name)
    {standardGeneric('getParam')})

#' @export
setGeneric('setParam', function(object, name, value)
    {standardGeneric('setParam')})

##################### Methods ####################

#' @importFrom methods slot
setMethod("getParam", "GapsParams", function(object, name)
{
    slot(object, name)
})

#' @importFrom methods slot<- validObject
setMethod("setParam", "GapsParams", function(object, name, value)
{
    slot(object, name) <- value
    validObject(object)
    return(object)
})

setMethod("show", "GapsParams", function(object)
{
    cat("CoGAPS Parameters Object\n")
    cat("Number of patterns to estimate: ", object@nFactor)
})

checkParamsWithData <- function(params, nRD, nCD, nRS, nCS)
{
    if (nRD != nRS | nCD != nCS)
        stop('D and S matrices have mismached dimensions')
    if (params@whichMatrixFixed == 'A' & nrow(params@fixedPatterns) != nRD)
        stop('invalid fixed pattern length')
    if (params@whichMatrixFixed == 'P' & ncol(params@fixedPatterns) != nCD)
        stop('invalid fixed pattern length')
}
