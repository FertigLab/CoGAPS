setMethod("show", signature("CogapsParams"),
function(object)
{
    cat("An Object of class \"CogapsParams\"\n")
    cat("nPatterns          ", object@nPatterns, "\n")
    cat("nIterations        ", object@nIterations, "\n")
    cat("outputFrequency    ", object@outputFrequency, "\n")
    cat("nCores             ", object@nCores, "\n")
    cat("singleCell         ", object@singleCell, "\n")
    cat("seed               ", object@seed, "\n")
    cat("messages           ", object@messages, "\n")
    cat("checkpointInterval ", object@checkpointInterval, "\n")
    cat("checkpointOutFile  ", object@checkpointOutFile, "\n")
})

#' @rdname setParam-methods
#' @aliases setParam
setMethod("setParam", signature(object="CogapsParams"),
function(object, whichParam, value)
{
    slot(object, whichParam) <- value
    validObject(object)
    return(object)
})

#' @rdname getParam-methods
#' @aliases getParam
setMethod("getParam", signature(object="CogapsParams"),
function(object, whichParam)
{
    slot(object, whichParam)
})