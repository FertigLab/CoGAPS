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
#' @importClassesFrom S4Vectors character_OR_NULL
setClass("CogapsParams", slots = c(
    nPatterns = "numeric",
    nIterations = "numeric",
    alphaA = "numeric",
    alphaP = "numeric",
    maxGibbsMassA = "numeric",
    maxGibbsMassP = "numeric",
    seed = "numeric",
    singleCell = "logical",
    distributed = "character_OR_NULL",
    nSets = "numeric",
    cut = "numeric",
    minNS = "numeric",
    maxNS = "numeric"
))

#' constructor for CogapsParams
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
        .Object@distributed <- NULL
        .Object@cut <- .Object@nPatterns
        .Object@nSets <- 4
        .Object@minNS <- ceiling(.Object@nSets / 2)
        .Object@maxNS <- .Object@minNS + .Object@nSets

        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

## defines a valid parameters object
setValidity("CogapsParams",
    function(object)
    {
        if (object@nPatterns <= 0 | object@nPatterns %% 1 != 0)
            "number of patterns must be an integer greater than zero"
        if (object@nIterations <= 0 | object@nIterations %% 1 != 0)
            "number of iterations must be an integer greater than zero"
        if (object@alphaA  <= 0 | object@alphaP <= 0)
            "alpha parameter must be greater than zero"
        if (object@maxGibbsMassA  <= 0 | object@maxGibbsMassP <= 0)
            "maxGibbsMass must be greater than zero"
        if (object@seed <= 0 | object@seed %% 1 != 0)
            "random seed must be an integer greater than zero"
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

#' set the value of parameters for distributed CoGAPS
#' @export
#' @docType methods
#' @rdname setDistributedParams-methods
#'
#' @description these parameters  are interrelated so they must be set together
#' @param object an object of type CogapsParams
#' @param cut a distributed CoGAPS parameter 
#' @param minNS a distributed CoGAPS parameter 
#' @param maxNS a distributed CoGAPS parameter 
#' @return the modified params object
#' @examples
#'  params <- new("CogapsParams")
#'  params <- setDistributedParams(3, 2, 4)
setGeneric("setDistributedParams", function(object, cut=NULL,
minNS=NULL, maxNS=NULL)
    {standardGeneric("setDistributedParams")})

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