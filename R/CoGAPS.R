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
CoGAPS <- function(D, S, nFactor=7, nEquil=1000, nSample=1000, nOutputs=1000,
nSnapshots=0, alphaA=0.01, alphaP=0.01, maxGibbmassA=100, maxGibbmassP=100,
seed=-1, messages=TRUE, singleCellRNASeq=FALSE, whichMatrixFixed = 'N',
fixedPatterns = matrix(0), checkpointInterval=0, ...)
{
    # get v2 arguments
    oldArgs <- list(...)
    if (!is.null(oldArgs$nOutR))
        nOutputs <- oldArgs$nOutR
    if (!is.null(oldArgs$max_gibbmass_paraA))
        maxGibbmassA <- oldArgs$max_gibbmass_paraA
    if (!is.null(oldArgs$max_gibbmass_paraP))
        maxGibbmassP <- oldArgs$max_gibbmass_paraP
    if (!is.null(oldArgs$sampleSnapshots) & is.null(oldArgs$numSnapshots))
        nSnapshots <- 100
    if (!is.null(oldArgs$sampleSnapshots) & !is.null(oldArgs$numSnapshots))
        nSnapshots <- oldArgs$numSnapshots
    if (missing(D) & !is.null(oldArgs$data))
        D <- oldArgs$data
    if (missing(S) & !is.null(oldArgs$unc))
        S <- oldArgs$unc

    # check arguments
    if (class(D) != "matrix" | class(S) != "matrix")
        stop('D and S must be matrices')
    if (any(D < 0) | any(S < 0))
        stop('D and S matrix must be non-negative')

    # run algorithm with call to C++ code
    nFactor <- floor(nFactor)
    nEquil <- floor(nEquil)
    nSample <- floor(nSample)
    result <- cogaps_cpp(D, S, nFactor, nEquil, floor(nEquil/10), nSample, nOutputs, nSnapshots,
        alphaA, alphaP, maxGibbmassA, maxGibbmassP, seed, messages,
        singleCellRNASeq, whichMatrixFixed, fixedPatterns, checkpointInterval)

    # backwards compatible with v2
    return(v2CoGAPS(result, ...))
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
#' @return list with A and P matrix estimates
#' @importFrom methods new
#' @export
gapsRun <- function(D, S, ABins=data.frame(), PBins=data.frame(), nFactor=7,
simulation_id="simulation", nEquil=1000, nSample=1000, nOutR=1000,
output_atomic=FALSE, fixedBinProbs=FALSE, fixedDomain="N", sampleSnapshots=TRUE,
numSnapshots=100, alphaA=0.01, nMaxA=100000, max_gibbmass_paraA=100.0,
alphaP=0.01, nMaxP=100000, max_gibbmass_paraP=100.0, seed=-1, messages=TRUE)
{
    #warning('gapsRun is deprecated with v3.0, use CoGAPS')
    CoGAPS(D, S, nFactor=nFactor, nEquil=nEquil, nSample=nSample, nOutputs=nOutR,
        nSnapshots=ifelse(sampleSnapshots,numSnapshots,0), alphaA=alphaA,
        alphaP=alphaP, maxGibbmassA=max_gibbmass_paraA, messages=messages,
        maxGibbmassP=max_gibbmass_paraP, seed=seed)
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
gapsMapRun <- function(D, S, FP, ABins=data.frame(), PBins=data.frame(),
nFactor=5, simulation_id="simulation", nEquil=1000, nSample=1000, nOutR=1000,
output_atomic=FALSE, fixedMatrix="P", fixedBinProbs=FALSE, fixedDomain="N",
sampleSnapshots=TRUE, numSnapshots=100, alphaA=0.01, nMaxA=100000,
max_gibbmass_paraA=100.0, alphaP=0.01, nMaxP=100000, max_gibbmass_paraP=100.0,
seed=-1, messages=TRUE)
{
    #warning('gapsMapRun is deprecated with v3.0, use CoGaps')
    CoGAPS(D, S, nFactor=nFactor, nEquil=nEquil, nSample=nSample, nOutputs=nOutR,
        nSnapshots=ifelse(sampleSnapshots,numSnapshots,0), alphaA=alphaA,
        alphaP=alphaP, maxGibbmassA=max_gibbmass_paraA, messages=messages,
        maxGibbmassP=max_gibbmass_paraP, seed=seed, whichMatrixFixed='P',
        fixedPatterns=as.matrix(FP))
}

v2CoGAPS <- function(result, ...)
{
    if (!is.null(list(...)$GStoGenes))
    {

    }
    return(result)
}