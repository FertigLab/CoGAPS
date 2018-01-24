#' CoGAPS Matrix Factorization Algorithm
#' 
#' @details calls the C++ MCMC code and performs Bayesian
#'matrix factorization returning the two matrices that reconstruct
#'the data matrix
#' @param D data matrix
#' @param S uncertainty matrix (std devs for chi-squared of Log Likelihood)
#' @return list with A and P matrix estimates
#' @importFrom methods new
#' @export
CoGAPS <- function(D, S, ...)
{
    # process v2 style arguments for backwards compatibility
    extraArgs <- list(...)
    params <- oldParams(new('GapsParams', 7, 1000, 1000), extraArgs)

    # check data
    if (missing(D) & length(extraArgs$data))
        D <- extraArgs$data
    if (missing(S) & length(extraArgs$unc))
        S <- extraArgs$unc
    if (any(D < 0) | any(S < 0))
        stop('D and S matrix must be non-negative')
    checkParamsWithData(params, nrow(D), ncol(D), nrow(S), ncol(S))

    # run algorithm with call to C++ code
    result <- cogaps_cpp(as.matrix(D), as.matrix(S),
        as.matrix(params$fixedPatterns), params)

    # backwards compatible with v2
    if (!is.null(extraArgs$GStoGenes))
    {
        #warning('the GStoGenes argument is deprecated with v3.0')
        if (is.null(extraArgs$plot) | extraArgs$plot)
            plotGAPS(result$Amean, result$Pmean)
        if (is.null(extraArgs$nPerm))
            extraArgs$nPerm <- 500
        GSP <- calcCoGAPSStat(result$Amean, result$Asd, extraArgs$GStoGenes,
            extraArgs$nPerm)
        result <- list(meanChi2=result$meanChi2, D=D, Sigma=S,
            Amean=result$Amean, Asd=result$Asd, Pmean=result$Pmean,
            Psd=result$Psd, GSUpreg=GSP$GSUpreg, GSDownreg=GSP$GSDownreg,
            GSActEst=GSP$GSActEst)
    }
    return(result)
}

#' Restart CoGAPS from Checkpoint File
#'
#' @details loads the state of a previous CoGAPS run from a file and
#'  continues the run from that point
#' @param D data matrix
#' @param S uncertainty matrix
#' @param path path to checkpoint file
#' @return list with A and P matrix estimates
#' @export
CoGapsFromCheckpoint <- function(D, S, path)
{
    cogapsFromCheckpoint_cpp(D, S, path)
}

GWCoGapsFromCheckpoint <- function(fname)
{
    #TODO
}

#' Display Information About Package Compilation
#' @export
displayBuildReport <- function()
{
    displayBuildReport_cpp()
}

#' Backwards Compatibility with v2
#'
#' @param D data matrix
#' @param S uncertainty matrix
#' @param ... v2 style parameters
#' @return list with A and P matrix estimates
#' @importFrom methods new
#' @export
gapsRun <- function(D, S, ...)
{
    #warning('gapsRun is deprecated with v3.0, use CoGAPS')
    params <- oldParams(new('GapsParams', 7, 1000, 1000), list(...))
    CoGAPS(D, S, params)
}

#' Backwards Compatibility with v2
#'
#' @param D data matrix
#' @param S uncertainty matrix
#' @param FP data.frame with rows giving fixed patterns for P
#' @param ... v2 style parameters
#' @return list with A and P matrix estimates
#' @importFrom methods new
#' @export
gapsMapRun <- function(D, S, FP, ...)
{
    #warning('gapsMapRun is deprecated with v3.0, use CoGaps')
    params <- oldParams(new('GapsParams', 7, 1000, 1000), list(...))
    CoGAPS(D, S, params)
}

# helper function to process v2 parameters
oldParams <- function(params, args)
{
    # standard params
    if (length(args$nFactor))      params@nFactor    <- as.integer(args$nFactor)
    if (length(args$nEquil))       params@nEquil     <- as.integer(args$nEquil)
    if (length(args$nSample))      params@nSample    <- as.integer(args$nSample)
    if (length(args$alphaA))       params@alphaA     <- args$alphaA
    if (length(args$alphaP))       params@alphaP     <- args$alphaP
    if (length(args$seed))         params@seed       <- as.integer(args$seed)
    if (length(args$messages))     params@messages   <- args$messages
    if (length(args$numSnapshots))
        params@nSnapshots <- as.integer(args$numSnapshots)

    # gapsMap params
    if (length(args$FP))
    {
        params@fixedPatterns <- as.matrix(args$FP)
        params@whichMatrixFixed <- 'P'
    }

    # return v3 style GapsParams
    return(params)    
}