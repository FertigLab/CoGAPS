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
#' @slot sparseOptimization speeds up performance with sparse data
#' (roughly >80% of data is zero), note this can only be used with the
#' default uncertainty
#' @slot distributed either "genome-wide" or "single-cell" indicating which
#' distributed algorithm should be used
#' @slot nSets [distributed parameter] number of sets to break data into
#' @slot cut [distributed parameter] number of branches at which to cut
#' dendrogram used in pattern matching
#' @slot minNS [distributed parameter] minimum of individual set contributions
#' a cluster must contain
#' @slot maxNS [distributed parameter] maximum of individual set contributions
#' a cluster can contain
#' @slot explicitSets [distributed parameter] specify subsets by index or name
#' @slot samplingAnnotation [distributed parameter] specify categories along
#' the rows (cols) to use for weighted sampling
#' @slot samplingWeight [distributed parameter] weights associated with 
#' samplingAnnotation
#' @slot subsetIndices set of indices to use from the data
#' @slot subsetDim which dimension (1=rows, 2=cols) to subset
#' @slot geneNames vector of names of genes in data
#' @slot sampleNames vector of names of samples in data
#' @slot fixedPatterns fix either 'A' or 'P' matrix to these values, in the
#' context of distributed CoGAPS (GWCoGAPS/scCoGAPS), the first phase is
#' skipped and fixedPatterns is used for all sets - allowing manual pattern
#' matching, as well as fixed runs of standard CoGAPS
#' @slot whichMatrixFixed either 'A' or 'P', indicating which matrix is fixed
#' @slot takePumpSamples whether or not to take PUMP samples
#' @slot checkpointInterval how many iterations between each checkpoint (set to 0 to disable)
#' @slot checkpointInFile file path to load checkpoint from
#' @slot checkpointOutFile file path where checkpoint should be written to
#' @importClassesFrom S4Vectors character_OR_NULL
setClass("CogapsParams", slots = c(
    nPatterns = "numeric",
    nIterations = "numeric",
    alphaA = "numeric",
    alphaP = "numeric",
    maxGibbsMassA = "numeric",
    maxGibbsMassP = "numeric",
    seed = "numeric",
    sparseOptimization = "logical",
    distributed = "character_OR_NULL",
    nSets = "numeric",
    cut = "numeric",
    minNS = "numeric",
    maxNS = "numeric",
    explicitSets = "ANY",
    samplingAnnotation = "character_OR_NULL",
    samplingWeight = "ANY",
    subsetIndices="ANY",
    subsetDim="numeric",
    geneNames="character_OR_NULL",
    sampleNames="character_OR_NULL",
    checkpointInterval="numeric",
    checkpointInFile="character_OR_NULL",
    checkpointOutFile="character_OR_NULL",
    fixedPatterns="ANY",
    whichMatrixFixed="character",
    takePumpSamples="logical"
))

#' constructor for CogapsParams
#' @param .Object CogapsParams object
#' @param distributed either "genome-wide" or "single-cell" indicating which
#' distributed algorithm should be used
#' @param ... initial values for slots
#' @return initialized CogapsParams object
#' @importFrom methods callNextMethod
setMethod("initialize", "CogapsParams",
    function(.Object, distributed=NULL, ...)
    {
        getMilliseconds <- function(time) floor((time$sec %% 1) * 1000)

        if (!is.null(list(...)$nSets))
            stop("nSets must be set after CogapsParams are intialized")
        if (!is.null(list(...)$cut))
            stop("cut must be set after CogapsParams are intialized")
        if (!is.null(list(...)$minNS))
            stop("minNS must be set after CogapsParams are intialized")
        if (!is.null(list(...)$maxNS))
            stop("maxNS must be set after CogapsParams are intialized")
        if (!is.null(distributed))
            if (distributed == "none") # allows it to be a pure string parameter
                distributed <- NULL
        .Object@distributed <- distributed
        
        .Object@nPatterns <- 7
        .Object@nIterations <- 50000
        .Object@alphaA <- 0.01
        .Object@alphaP <- 0.01
        .Object@maxGibbsMassA <- 100
        .Object@maxGibbsMassP <- 100
        .Object@seed <- getMilliseconds(as.POSIXlt(Sys.time()))
        .Object@sparseOptimization <- FALSE
        .Object@cut <- .Object@nPatterns
        .Object@nSets <- 4
        .Object@minNS <- ceiling(.Object@nSets / 2)
        .Object@maxNS <- .Object@minNS + .Object@nSets
        .Object@explicitSets <- NULL
        .Object@samplingAnnotation <- NULL
        .Object@samplingWeight <- NULL
        .Object@subsetIndices <- NULL
        .Object@subsetDim <- 0
        .Object@geneNames <- NULL
        .Object@sampleNames <- NULL
        .Object@fixedPatterns <- NULL
        .Object@whichMatrixFixed <- 'N'
        .Object@takePumpSamples <- FALSE
        .Object@checkpointInterval <- 0
        .Object@checkpointInFile <- NULL
        .Object@checkpointOutFile <- NULL

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
        if (object@alphaA <= 0 | object@alphaP <= 0)
            "alpha parameter must be greater than zero"
        if (object@maxGibbsMassA <= 0 | object@maxGibbsMassP <= 0)
            "maxGibbsMass must be greater than zero"
        if (object@seed <= 0 | object@seed %% 1 != 0)
            "random seed must be an integer greater than zero"

        if (!(object@whichMatrixFixed %in% c("A", "P", "N")))
            stop("Invalid choice of fixed matrix, must be 'A' or 'P'")
        if (!is.null(object@fixedPatterns) & object@whichMatrixFixed == "N")
            stop("fixedPatterns passed without setting whichMatrixFixed")
        if (object@whichMatrixFixed %in% c("A", "P") & is.null(object@fixedPatterns))
            stop("whichMatrixFixed is set without passing fixedPatterns")

        if (!(object@subsetDim %in% c(0,1,2)))
            stop("invalid subset dimension")
        if (object@subsetDim > 0 & is.null(object@subsetIndices))
            stop("subsetDim provided without subsetIndices")

        if (!is.null(object@distributed))
        {
            if (!(object@distributed %in% c("genome-wide", "single-cell")))
                "distributed method must be either 'genome-wide' or 'single-cell'"
            if (!is.null(object@fixedPatterns) & is.null(object@explicitSets))
                "doing manual pattern matching without using explicit subsets"
            if (object@distributed == "single-cell" & object@whichMatrixFixed == "P")
                "can't fix P matrix when running single-cell CoGAPS"
            if (object@distributed == "genome-wide" & object@whichMatrixFixed == "A")
                "can't fix A matrix when running genome-wide CoGAPS"
            if (object@minNS <= 1 | object@minNS %% 1 != 0)
                "minNS must be an integer greater than one"
            if (object@nSets <= 1 | object@nSets %% 1 != 0)
                "minNS must be an integer greater than one"
            if (length(unique(object@samplingAnnotation)) != length(object@samplingWeight))
                "samplingWeight has mismatched size with amount of distinct annotations"
            if (object@cut > object@nPatterns)
                "cut must be less than or equal to nPatterns"
            if (length(object@samplingWeight) & is.null(names(object@samplingWeight)))
                "samplingWeight must be a named vector"

            if (!is.null(object@explicitSets))
            {
                if (!is(object@explicitSets, "list"))
                    "explicitSets must be a list"
                if (length(object@explicitSets) != object@nSets)
                    "nSets doesn't match length of explicitSets"
                if (!is.null(object@samplingAnnotation))
                    "explicitSets and samplingAnnotation/samplingWeight are both set"
                isNum <- sapply(object@explicitSets, function(s) is(s, "numeric"))
                isChar <- sapply(object@explicitSets, function(s) is(s, "character"))
                if (!all(isNum) & !all(isChar))
                    "explicitSets must be a list of numeric or character"
            }
        }
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
#' @param nSets number of sets to break data into
#' @param cut number of branches at which to cut dendrogram used in
#' pattern matching
#' @param minNS minimum of individual set contributions a cluster must contain
#' @param maxNS maximum of individual set contributions a cluster can contain
#' @return the modified params object
#' @examples
#'  params <- new("CogapsParams")
#'  params <- setDistributedParams(params, 5)
setGeneric("setDistributedParams", function(object, nSets=NULL, cut=NULL,
minNS=NULL, maxNS=NULL)
    {standardGeneric("setDistributedParams")})

#' set the annotation labels and weights for subsetting the data
#' @export
#' @docType methods
#' @rdname setAnnotationWeights-methods
#'
#' @description these parameters are interrelated so they must be set together
#' @param object an object of type CogapsParams
#' @param annotation vector of labels
#' @param weights vector of weights
#' @return the modified params object
#' @examples
#'  params <- new("CogapsParams")
#'  params <- setAnnotationWeights(params, c('a', 'b', 'c'), c(1,2,1))
setGeneric("setAnnotationWeights", function(object, annotation, weights)
    {standardGeneric("setAnnotationWeights")})

#' set the fixed patterns for either the A or the P matrix
#' @export
#' @docType methods
#' @rdname setFixedPatterns-methods
#'
#' @description these parameters are interrelated so they must be set together
#' @param object an object of type CogapsParams
#' @param fixedPatterns values for either the A or P matrix
#' @param whichMatrixFixed either 'A' or 'P' indicating which matrix is fixed
#' @return the modified params object
#' @examples
#' params <- new("CogapsParams")
#' data(GIST)
#' params <- setFixedPatterns(params, getSampleFactors(GIST.result), 'P')
setGeneric("setFixedPatterns", function(object, fixedPatterns, whichMatrixFixed)
    {standardGeneric("setFixedPatterns")})

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
