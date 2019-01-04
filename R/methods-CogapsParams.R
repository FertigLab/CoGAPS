setMethod("show", signature("CogapsParams"),
function(object)
{
    cat("-- Standard Parameters --\n")
    cat("nPatterns           ", object@nPatterns, "\n")
    cat("nIterations         ", object@nIterations, "\n")
    cat("seed                ", object@seed, "\n")
    cat("singleCell          ", object@singleCell, "\n")
    cat("sparseOptimization  ", object@sparseOptimization, "\n")
    if (!is.null(object@distributed))
        cat("distributed         ", object@distributed, "\n")
    cat("\n")
    cat("-- Sparsity Parameters --\n")
    if (object@alphaA == object@alphaP)
    {
        cat("alpha         ", object@alphaA, "\n")
    }
    else
    {
        cat("alphaA        ", object@alphaA, "\n")
        cat("alphaP        ", object@alphaP, "\n")
    }

    if (object@maxGibbsMassA == object@maxGibbsMassP)
    {
        cat("maxGibbsMass  ", object@maxGibbsMassA, "\n")
    }
    else
    {
        cat("maxGibbsMassA ", object@maxGibbsMassA, "\n")
        cat("maxGibbsMassP ", object@maxGibbsMassP, "\n")
    }
    if (!is.null(object@distributed))
    {
        cat("\n")
        cat("-- Distributed CoGAPS Parameters --", "\n")
        cat("nSets         ", object@nSets, "\n")
        cat("cut           ", object@cut, "\n")
        cat("minNS         ", object@minNS, "\n")
        cat("maxNS         ", object@maxNS, "\n")
    }
})

#' @rdname setParam-methods
#' @aliases setParam
#' @importFrom methods slot slot<- validObject
setMethod("setParam", signature(object="CogapsParams"),
function(object, whichParam, value)
{
    if (whichParam == "alpha")
    {
        object@alphaA <- value
        object@alphaP <- value
    }
    else if (whichParam == "maxGibbsMass")
    {
        object@maxGibbsMassA <- value
        object@maxGibbsMassP <- value
    }
    else if (whichParam %in% c("nSets", "cut", "minNS", "maxNS"))
    {
        stop("please set this parameter with setDistributedParams")
    }
    else if (whichParam %in% c("samplingAnnotation", "samplingWeight"))
    {
        stop("please set this parameter with setAnnotationWeights")
    }
    else if (whichParam == "nPatterns")
    {
        object@nPatterns <- value
        object@cut <- min(object@cut, object@nPatterns)
    }
    else
    {
        slot(object, whichParam) <- value
    }
    validObject(object)
    return(object)
})

#' @rdname setDistributedParams-methods
#' @aliases setDistributedParams
#' @importFrom methods slot
setMethod("setDistributedParams", signature(object="CogapsParams"),
function(object, nSets, cut, minNS, maxNS)
{
    message("setting distributed parameters - call this again if you change ",
        "nPatterns")

    object@nSets <- nSets

    object@cut <- ifelse(is.null(cut), object@nPatterns, cut)
    object@minNS <- ifelse(is.null(minNS), ceiling(object@nSets / 2), minNS)
    object@maxNS <- ifelse(is.null(maxNS), object@minNS + object@nSets, maxNS)

    validObject(object)
    return(object)
})

#' @rdname setAnnotationWeights-methods
#' @aliases setAnnotationWeights
setMethod("setAnnotationWeights", signature(object="CogapsParams"),
function(object, annotation, weights)
{
    object@samplingAnnotation <- annotation
    object@samplingWeight <- weights

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