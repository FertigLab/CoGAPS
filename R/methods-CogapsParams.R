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

#' @rdname parseOldParams-methods
#' @aliases parseOldParams
setMethod("parseOldParams", signature(object="CogapsParams"),
function(object, oldArgs)
{
    helper <- function(arg, params, newArg)
    {
        if (!is.null(oldArgs[[arg]]))
        {
            warning(paste("parameter", arg, "is deprecated, it will still",
                "work, but setting", newArg, "in the params object is the",
                "preferred method"))
            params <- setParam(params, newArg, oldArgs[[arg]])
            oldArgs[[arg]] <- NULL
        }            
        return(params)
    }

    object <- helper("nFactor", object, "nPatterns")
    object <- helper("nIter", object, "nIterations")
    object <- helper("nEquil", object, "nIterations")
    object <- helper("nSample", object, "nIterations")
    object <- helper("nOutR", object, "outputFrequency")
    object <- helper("nOutput", object, "outputFrequency")
    object <- helper("maxGibbmassA", object, "maxGibbsMassA")
    object <- helper("max_gibbmass_paraA", object, "maxGibbsMassA")
    object <- helper("maxGibbmassP", object, "maxGibbsMassP")
    object <- helper("max_gibbmass_paraP", object, "maxGibbsMassP")
    object <- helper("singleCellRNASeq", object, "singleCell")

    if (!is.null(oldArgs$nSnapshots) | !is.null(oldArgs$sampleSnapshots) | !is.null(oldArgs$numSnapshots))
    {
        warning("snapshots not currently supported in release build")
        oldArgs$nSnapshots <- NULL
        oldArgs$sampleSnapshots <- NULL
        oldArgs$numSnapshots <- NULL
    }
    if (!is.null(oldArgs$fixedPatterns))
        stop("pass fixed matrix in with 'fixedMatrix' argument")
    if (!is.null(oldArgs$S))
        stop("pass uncertainty matrix in with 'uncertainty', not 'S'")

    return(object)
})

#' @rdname parseDirectParams-methods
#' @aliases parseDirectParams
setMethod("parseDirectParams", signature(object="CogapsParams"),
function(object, args)
{
    for (s in slotNames(object))
    {
        if (!is.null(args[[s]]))
        {
            object <- setParam(object, s, args[[s]])
        }
    }
    return(object)
})