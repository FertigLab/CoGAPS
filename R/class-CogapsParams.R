#' CogapsParams
#' @export 
#' @rdname CogapsParams-class
#'
#' @description Encapsulates all parameters for the CoGAPS algorithm
#' @slot nPatterns number of patterns CoGAPS will learn
#' @slot nIterations number of iterations for each phase of the algorithm
#' @slot alphaA sparsity parameter for feature matrix
#' @slot alphaP sparsity parameter for sample matrix
#' @slot maxGibbsMassA atomic mass restriction for feature matrix
#' @slot maxGibbsMassP atomic mass restriction for sample matrix
#' @slot seed random number generator seed
#' @slot singleCell is the data single cell?
#' @slot distributed either "genome-wide" or "single-cell" indicating which
#' distributed algorithm should be used
#' @slot nSets [distributed parameter] number of sets to break data into
#' @slot cut [distributed parameter] number of branches at which to cut
#' dendrogram used in pattern matching
#' @slot minNS [distributed parameter] minimum of individual set contributions
#' a cluster must contain
#' @slot maxNS [distributed parameter] maximum of individual set contributions
#' a cluster can contain
setClass("CogapsParams", slots = c(
    nPatterns = "numeric",
    nIterations = "numeric",
    alphaA = "numeric",
    alphaP = "numeric",
    maxGibbsMassA = "numeric",
    maxGibbsMassP = "numeric",
    seed = "numeric",
    singleCell = "logical",
    distributed = "character",
    cut = "numeric",
    nSets = "numeric",
    minNS = "numeric",
    maxNS = "numeric"
))

#' Constructor for CogapsParams
#' @return initialized CogapsParams object
#' @importFrom methods callNextMethod
setMethod("initialize", "CogapsParams",
    function(.Object, ...)
    {
        getMilliseconds <- function(time) floor((time$sec %% 1) * 1000)

        .Object@nPatterns <- 7
        .Object@nIterations <- 1000
        .Object@alphaA <- 0.01
        .Object@alphaP <- 0.01
        .Object@maxGibbsMassA <- 100
        .Object@maxGibbsMassP <- 100
        .Object@seed <- getMilliseconds(as.POSIXlt(Sys.time()))
        .Object@singleCell <- FALSE
        .Object@distributed <- ""
        .Object@cut <- .Object@nPatterns
        .Object@nSets <- 3
        .Object@minNS <- ceiling(.Object@nSets / 2)
        .Object@maxNS <- .Object@minNS + .Object@nSets

        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

setValidity("CogapsParams",
    function(object)
    {
        if (object@nPatterns <= 0 | object@nPatterns %% 1 != 0)
            "number of patterns must be an integer greater than zero"
        if (object@nIterations <= 0 | object@nIterations %% 1 != 0)
            "number of iterations must be an integer greater than zero"
        if (object@outputFrequency <= 0 | object@outputFrequency %% 1 != 0)
            "the output frequency must be an integer greater than zero"
        if (object@alphaA  <= 0 | object@alphaP <= 0)
            "alpha parameter must be greater than zero"
        if (object@maxGibbsMassA  <= 0 | object@maxGibbsMassP <= 0)
            "maxGibbsMass must be greater than zero"
        if (object@seed <= 0 | object@seed %% 1 != 0)
            "random seed must be an integer greater than zero"
        if (!(object@whichMatrixFixed %in% c("N", "A", "P")))
            "whichMatrixFixed must be either A or P (N in the case of neither)"
        if (object@checkpointInterval <= 0 | object@checkpointInterval %% 1 != 0)
            "checkpointInterval must be an integer greater than zero"
        if (object@nCores <= 0 | object@nCores %% 1 != 0)
            "number of cores must be an integer greater than zero"
        if (object@minNS <= 1 | object@minNS %% 1 != 0)
            "minNS must be an integer greater than one"
    }
)

#' set the value of a parameter
#' @export
#' @docType methods
#' @rdname setParam-methods
#'
#' @param object an object of type CogapsParams
#' @param whichParam a string with the name of the parameter to be changed
#' @param value the value to set the parameter to
#' @return the modified params object
#' @examples
#'  params <- new("CogapsParams")
#'  params <- setParam(params, "seed", 123)
setGeneric("setParam", function(object, whichParam, value)
    {standardGeneric("setParam")})

#' get the value of a parameter
#' @export
#' @docType methods
#' @rdname getParam-methods
#'
#' @param object an object of type CogapsParams
#' @param whichParam a string with the name of the requested parameter
#' @return the value of the parameter
#' @examples
#'  params <- new("CogapsParams")
#'  getParam(params, "seed")
setGeneric("getParam", function(object, whichParam)
    {standardGeneric("getParam")})

#' parse list of old-style parameters, store relevant values
#' @docType methods
#' @rdname parseOldParams-methods
#'
#' @param object an object of type CogapsParams
#' @param oldArgs named list of deprecated arguments
#' @return an object of type CogapsParams
setGeneric("parseOldParams", function(object, oldArgs)
    {standardGeneric("parseOldParams")})

#' parse list of parameters passed directly to CoGAPS
#' @docType methods
#' @rdname parseDirectParams-methods
#'
#' @param object an object of type CogapsParams
#' @param oldArgs named list of arguments
#' @return an object of type CogapsParams
setGeneric("parseDirectParams", function(object, args)
    {standardGeneric("parseDirectParams")})

