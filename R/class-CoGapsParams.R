#' CoGapsParams
#' @export 
#'
#' @description Encapsulates all parameters for the CoGAPS algorithm
setClass("CoGapsParams", slots = c(
    nFactor = "numeric",
    nIter = "numeric",
    outputFrequency = "numeric",
    alphaA = "numeric",
    alphaP = "numeric",
    maxGibbsMassA = "numeric",
    maxGibbsMassP = "numeric",
    seed = "numeric",
    messages = "logical",
    singleCellRNASeq = "logical",
    whichMatrixFixed = "character",
    checkpointInterval = "numeric",
    checkpointOutFile = "character", 
    nCores = "numeric"
))

#' Constructor for CoGapsParams
#' @return initialized CoGapsParams object
#' @importFrom methods callNextMethod
setMethod("initialize", "CoGapsParams",
    function(.Object, ...)
    {
        getMilliseconds <- function(time) floor((time$sec %% 1) * 1000)

        .Object@nFactor <- 7
        .Object@nIter <- 1000
        .Object@outputFrequency <- 500
        .Object@alphaA <- 0.01
        .Object@alphaP <- 0.01
        .Object@maxGibbsMassA <- 100
        .Object@maxGibbsMassP <- 100
        .Object@seed <- getMilliseconds(as.POSIXlt(Sys.time()))
        .Object@messages <- TRUE
        .Object@singleCellRNASeq <- FALSE
        .Object@whichMatrixFixed <- "N"
        .Object@checkpointInterval <- 0
        .Object@checkpointOutFile <- "gaps_checkpoint.out"
        .Object@nCores <- 1

        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

setValidity("CoGapsParams",
    function(object)
    {
        if (object@nFactor <= 0 || object@nFactor %% 1 != 0)
            "number of patterns must be an integer greater than zero"
        if (object@nIter <= 0 || object@nIter %% 1 != 0)
            "number of iterations must be an integer greater than zero"
        if (object@outputFrequency <= 0 || object@outputFrequency %% 1 != 0)
            "the output frequency must be an integer greater than zero"
        if (object@alphaA  <= 0 || object@alphaP <= 0)
            "alpha parameter must be greater than zero"
        if (object@maxGibbsMassA  <= 0 || object@maxGibbsMassP <= 0)
            "maxGibbsMass must be greater than zero"
        if (object@seed <= 0 || object@seed %% 1 != 0)
            "random seed must be an integer greater than zero"
        if (!(object@whichMatrixFixed %in% c("N", "A", "P")))
            "whichMatrixFixed must be either A or P (N in the case of neither)"
        if (object@checkpointInterval <= 0 || object@checkpointInterval %% 1 != 0)
            "checkpointInterval must be an integer greater than zero"
        if (object@nCores <= 0 || object@nCores %% 1 != 0)
            "number of cores must be an integer greater than zero"
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
setGeneric("setParam", function(object, whichParam, value)
    {standardGeneric("setParam")})

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
setGeneric("getParam", function(object, whichParam)
    {standardGeneric("getParam")})

#' @rdname setParam-methods
#' @aliases setParam
setMethod("setParam", signature(object="CoGapsParams"), 
    function(object, whichParam, value)
    {
        slot(params, whichParam) <- value
        validObject(params)
        return(params)
    }
)

#' @rdname getParam-methods
#' @aliases getParam
setMethod("getParam", signature(object="CoGapsParams"), 
    function(object, whichParam)
    {
        slot(params, whichParam)
    }
)
