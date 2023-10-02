#' CogapsParams constructor
#' @export
#'
#' @description create a CogapsParams object
#' @param ... parameters for the initialization method
#' @return CogapsParams object
#' @examples
#' params <- CogapsParams(nPatterns=10)
#' params
CogapsParams <- function(...)
{
    if (!is.null(list(...)$singleCell))
        warning("singleCell has been deprecated, this parameter will be ignored")
    removeDeprecated <- function(..., singleCell) new("CogapsParams", ...)
    removeDeprecated(...)
}

setMethod("show", signature("CogapsParams"),
function(object)
{
    cat("-- Standard Parameters --\n")
    cat("nPatterns           ", object@nPatterns, "\n")
    cat("nIterations         ", object@nIterations, "\n")
    cat("seed                ", object@seed, "\n")
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

    if (!is.null(object@geneNames))
    {
        cat("\n")
        cat(length(object@geneNames), "gene names provided\n")
        cat("first gene name:", object@geneNames[1], "\n")
    }

    if (!is.null(object@sampleNames))
    {
        cat("\n")
        cat(length(object@sampleNames), "sample names provided\n")
        cat("first sample name:", object@sampleNames[1], "\n")
    }
    if (!is.null(object@checkpointInFile) || !is.null(object@checkpointOutFile))
    {
      cat("\n")
      cat("-- Checkpoint parameters--", "\n")
      cat("checkpointInterval          ", object@checkpointInterval, "\n")
      if(object@checkpointInterval == 0){
        cat("Warning!! Setting checkpointInterval=0 disables checkpoint logging.", "\n")
      }
      if(!is.null(object@checkpointInFile)){
        cat("checkpointInFile          ", checkpointInFile, "\n")
      }
      if (!is.null(object@checkpointOutFile)){
        cat("checkpointOutFile          ", object@checkpointOutFile, "\n")
      }
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
        stop("please set \'", whichParam, "\' with setDistributedParams")
    }
    else if (whichParam %in% c("samplingAnnotation", "samplingWeight"))
    {
        stop("please set \'", whichParam, "\' with setAnnotationWeights")
    }
    else if (whichParam %in% c("fixedPatterns", "whichMatrixFixed"))
    {
        stop("please set \'", whichParam, "\' with setFixedPatterns")
    }
    else if (whichParam == "nPatterns")
    {
        object@nPatterns <- value
        object@cut <- min(object@cut, object@nPatterns)
    }
    else if (whichParam == "checkpointInFile")
    {
      object@checkpointInFile <- value
    }
    else if (whichParam == "checkpointOutFile")
    {
      object@checkpointOutFile <- value
    }
    else if (whichParam == "checkpointInterval")
    {
      object@checkpointInterval <- value
      if(value==0){
        object@checkpointOutFile <- NULL 
        stop("setting checkpointInterval=0 disables checkpoint writing")
      }
    }
    else if (whichParam == "distributed")
    {
        if (value == "none")
            object@distributed <- NULL
        else
            object@distributed <- value
    }
    else if (whichParam %in% c("singleCell"))
    {
        warning(whichParam, " has been deprecated, this parameter will be ignored")
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

    object@nSets <- ifelse(is.null(nSets), object@nSets, nSets)
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

#' @rdname setFixedPatterns-methods
#' @aliases setFixedPatterns
setMethod("setFixedPatterns", signature(object="CogapsParams"),
function(object, fixedPatterns, whichMatrixFixed)
{
    object@fixedPatterns <- fixedPatterns
    object@whichMatrixFixed <- whichMatrixFixed

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