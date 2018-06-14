setClass("CoGapsParams", slots = c(
    nFactor = "numeric",
    nIter = "numeric",
    outputFrequency = "numeric",
    alphaA = "numeric",
    alphaP = "numeric",
    maxGibbmassA = "numeric",
    maxGibbMassP = "numeric",
    seed = "numeric",
    messages = "numeric",
    singleCellRNASeq = "logical",
    whichMatrixFixed = "character",
    fixedMatrix = "matrix",
    checkpointInterval = "numeric",
    checkpointFile = "character", 
    nCores = "numeric"
))

setMethod("initialize", "CoGapsParams",
    function(.Object, ...)
    {
        .Object@nFactor <- 7
        .Object@nIter <- 1000
        .Object@outputFrequency <- 1000
        .Object@alphaA <- 0.01
        .Object@alphaP <- 0.01
        .Object@maxGibbmassA <- 100
        .Object@maxGibbmassP <- 100
        .Object@seed <- -1
        .Object@messages <- TRUE
        .Object@singleCellRNASeq <- FALSE
        .Object@whichMatrixFixed <- 'N'
        .Object@fixedMatrix <- matrix(0)
        .Object@checkpointInterval <- 0
        .Object@checkpointFile <- "gaps_checkpoint.out"
        .Object@nCores <- 1

        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

setValidity("CoGapsParams",
    function(object)
    {
        if (object@nFactor < 0)
            stop('number of patterns must be non-negative')
        if (object@nEquil < 0)
            stop('number of iterations for burn-in must be non-negative')
        if (object@nSample < 0)
            stop('number of iterations for sampling must be non-negative')
        if (object@nOutputs < 0)
            stop('number of iterations to print status into R by iterations must be non-negative')
        if (object@nSnapshots < 0)
            stop('number of samples to capture must be non-negative')
        if (object@alphaA  < 0)
            stop('sparsity parameter for A domain must be non-negative')
        if (object@alphaP  < 0)
            stop('sparsity parameter for P domain must be non-negative')
        if (object@maxGibbmassA < 0)
            stop('limit must be non-negative')
        if (object@maxGibbmassP < 0)
            stop('limit must be non-negative')
        if (object@nCores < 0)
            stop('number of cpu cores to run in parallel order must be non-negative')
    }
)

#' set the value of a parameter
#' @export
#' @docType methods
#' @rdname setParam-methods
#'
#' @param params an object of type CoGapsParams
#' @param whichParam a string with the name of the parameter to be changed
#' @param value the value to set the parameter to
#' @return the modified params object
#' @examples
#'  params <- new("CoGapsParams")
#'  params <- setParam(params, "seed", 123)
setGeneric('setParam', function(params, whichParam, value)
    {standardGeneric('setParam')})

#' get the value of a parameter
#' @export
#' @docType methods
#' @rdname getParam-methods
#'
#' @param params an object of type CoGapsParams
#' @param whichParam a string with the name of the requested parameter
#' @return the value of the parameter
#' @examples
#'  params <- new("CoGapsParams")
#'  getParam(params, "seed")
setGeneric('getParam', function(params, whichParam)
    {standardGeneric('getParam')})

#' @rdname setParam-methods
#' @aliases setParam
setMethod('setParam', signature(params='CoGapsParams'), 
    function(params, whichParam, value)
    {
        if (whichParam == "seed")
            params@seed <- value
        if (whichParam == "nFactor")
            params@nFactor <- value

        # TODO(Hyejune) repeat this for all parameters

        return(params)
    }
)

#' @rdname getParam-methods
#' @aliases getParam
setMethod('getParam', signature(params='CoGapsParams'), 
    function(params, whichParam)
    {
        if (whichParam == "seed")
            return(params@seed)
        if (whichParam == "nFactor")
            return(params@nFactor)
        
        # TODO(Hyejune) repeat this for all parameters    
    }
)

