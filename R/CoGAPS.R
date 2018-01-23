#' CoGAPS Matrix Factorization Algorithm
#' 
#' @details calls the C++ MCMC code and performs Bayesian
#'matrix factorization returning the two matrices that reconstruct
#'the data matrix
#' @param D data matrix
#' @param S uncertainty matrix (std devs for chi-squared of Log Likelihood)
#' @param params GapsParams object 
#' @return list with A and P matrix estimates
#' @export
CoGaps <- function(D, S, params = new('GapsParams', 7, 1000, 1000), ...)
{
    # process v2 style arguments for backwards compatibility
    params <- oldParams(params, list(...))

    # check data
    checkParamsWithData(params, nrow(D), ncol(D), nrow(S), ncol(S))
    if (any(D) < 0 | any(S) < 0) # too slow for large matrices ?
        stop('D and S matrix must be non-negative')

    # run algorithm
    result <- cogaps_cpp(as.matrix(D), as.matrix(S), nFactor, alphaA,
        alphaP, nEquil, floor(nEquil/10), nSample, max_gibbmass_paraA,
        max_gibbmass_paraP, fixedPatterns, whichMatrixFixed, seed, messages,
        singleCellRNASeq, nOutR, numSnapshots, checkpoint_interval)

    # backwards compatible with v2
    if (length(list(...)$GStoGenes))
    {
        warning('the GStoGenes argument is depracted with v3.0,')
        if (list(...)$plot)
            plotGAPS(result$Amean, result$Pmean)
        GSP <- calcCoGAPSStat(result$Amean, result$Asd, list(...)$GStoGenes,
            list(...)$nPerm)
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
    # TODO add checksum to make sure D,S are correct
    cogapsFromCheckpoint_cpp(D, S, path)
}

#' Backwards Compatibility with v2
#'
#' @param D data matrix
#' @param S uncertainty matrix
#' @param ... v2 style parameters
#' @return list with A and P matrix estimates
#' @export
gapsRun <- function(D, S, ...)
{
    warning('gapsRun is depracated with v3.0, use CoGAPS')
    params <- new('GapsParams', 7, 1000, 1000)
    params <- oldParams(params, list(...))
    CoGAPS(D, S, params)
}

#' Backwards Compatibility with v2
#'
#' @param D data matrix
#' @param S uncertainty matrix
#' @param FP data.frame with rows giving fixed patterns for P
#' @param ... v2 style parameters
#' @return list with A and P matrix estimates
#' @export
gapsMapRun <- function(D, S, FP, ...)
{
    warning('gapsMapRun is depracated with v3.0, use CoGaps')
    params <- new('GapsParams', 7, 1000, 1000)
    params <- oldParams(params, list(...))
    CoGAPS(D, S, params)
}

# helper function to process v2 parameters
oldParams <- function(params, args)
{
    # standard params
    if (length(args$nFactor))      params@nFactor    <- args$nFactor
    if (length(args$nEquil))       params@nEquil     <- args$nEquil
    if (length(args$nSample))      params@nSample    <- args$nSample
    if (length(args$numSnapshots)) params@nSnapshots <- args$numSnapshots
    if (length(args$alphaA))       params@alphaA     <- args$alphaA
    if (length(args$alphaP))       params@alphaP     <- args$alphaP
    if (length(args$seed))         params@seed       <- args$seed
    if (length(args$messages))     params@messages   <- args$messages

    # gapsMap params
    if (length(args$FP))
    {
        params@fixedPatterns <- as.matrix(args$FP)
        params@whichMatrixFixed <- 'P'
    }

    # return v3 style GapsParams
    return(params)    
}
