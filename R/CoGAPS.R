# regenerate `SimpSim.result` with correct row/col names
# need to pass whole GSets list
# not one element
# smooth patterns nbd
# bluered not quotes - it's a function
# make sure import heatmap.2 from gplots
# make sure GWCoGAPS vignette is included
# look into conference poster/abstract submission
# write first draft of abstract

#' CoGAPS Matrix Factorization Algorithm
#' 
#' @details calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix
#' @param D data matrix
#' @param S uncertainty matrix (std devs for chi-squared of Log Likelihood)
#' @param nFactor number of patterns (basis vectors, metagenes), which must be
#' greater than or equal to the number of rows of FP
#' @param nEquil number of iterations for burn-in
#' @param nSample number of iterations for sampling
#' @param nOutputs how often to print status into R by iterations
#' @param nSnapshots the number of individual samples to capture
#' @param alphaA sparsity parameter for A domain
#' @param alphaP sparsity parameter for P domain
#' @param maxGibbmassA limit truncated normal to max size
#' @param maxGibbmassP limit truncated normal to max size
#' @param seed a positive seed is used as-is, while any negative seed tells
#' the algorithm to pick a seed based on the current time
#' @param messages display progress messages
#' @param singleCellRNASeq indicates if the data is single cell RNA-seq data
#' @param whichMatrixFixed character to indicate whether A or P matric contains
#' the fixed patterns
#' @param fixedPatterns matrix of fixed values in either A or P matrix
#' @param checkpointInterval time (in seconds) between creating a checkpoint
#' @param checkpointFile name of the checkpoint file
#' @param nCores number of cpu cores to run in parallel over
#' @param ... keeps backwards compatibility with arguments from older versions
#' @return list with A and P matrix estimates
#' @importFrom methods new
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D, SimpSim.S, nFactor=3, nOutputs=250)
#' @export
CoGAPS <- function(D, S, nFactor=7, nEquil=250, nSample=250, nOutputs=1000,
nSnapshots=0, alphaA=0.01, alphaP=0.01, maxGibbmassA=100, maxGibbmassP=100,
seed=NA, messages=TRUE, singleCellRNASeq=FALSE, whichMatrixFixed='N',
fixedPatterns=matrix(0), checkpointInterval=0, 
checkpointFile="gaps_checkpoint.out", nCores=1, ...)
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

    # get pump arguments - hidden for now from user
    pumpThreshold <- "unique"
    nPumpSamples <- 0
    if (!is.null(list(...)$pumpThreshold))
        pumpThreshold <- list(...)$pumpThreshold
    if (!is.null(list(...)$nPumpSamples))
        pumpThreshold <- list(...)$nPumpSamples

    # check arguments
    if (class(D) != "matrix" | class(S) != "matrix")
        stop('D and S must be matrices')
    if (any(D < 0) | any(S < 0))
        stop('D and S matrix must be non-negative')
    if (nrow(D) != nrow(S) | ncol(D) != ncol(S))
        stop('D and S matrix have different dimensions')
    if (whichMatrixFixed == 'A' & nrow(fixedPatterns) != nrow(D))
        stop('invalid number of rows for fixedPatterns')
    if (whichMatrixFixed == 'A' & ncol(fixedPatterns) > nFactor)
        stop('invalid number of columns for fixedPatterns')
    if (whichMatrixFixed == 'P' & nrow(fixedPatterns) > nFactor)
        stop('invalid number of rows for fixedPatterns')
    if (whichMatrixFixed == 'P' & ncol(fixedPatterns) != ncol(D))
        stop('invalid number of columns for fixedPatterns')
    thresholdEnum <- c("unique", "cut")
    
    # get seed
    if (is.na(seed))
    {
        seed <- 0 # TODO get time in milliseconds
    }

    # run algorithm with call to C++ code
    result <- cogaps_cpp(D, nFactor, nEquil, nOutputs, seed, alphaA, alphaP,
        maxGibbmassA, maxGibbmassP, messages, singleCellRNASeq, nCores)
    
    # label matrices and return list
    patternNames <- paste('Patt', 1:nFactor, sep='')
    rownames(result$Amean) <- rownames(result$Asd) <- rownames(D)
    colnames(result$Amean) <- colnames(result$Asd) <- patternNames
    rownames(result$Pmean) <- rownames(result$Psd) <- patternNames
    colnames(result$Pmean) <- colnames(result$Psd) <- colnames(D)
    return(v2CoGAPS(result, ...)) # backwards compatible with v2
}

#' Restart CoGAPS from Checkpoint File
#' @export
#'
#' @details loads the state of a previous CoGAPS run from a file and
#'  continues the run from that point
#' @param D data matrix
#' @param S uncertainty matrix
#' @param path path to checkpoint file
#' @param checkpointFile name for future checkpooints made
#' @return list with A and P matrix estimates
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D, SimpSim.S, nFactor=3, nOutputs=250)
CoGapsFromCheckpoint <- function(D, S, path, checkpointFile=NA)
{
    if (is.na(checkpointFile))
        checkpointFile <- path
    cogapsFromCheckpoint_cpp(D, S, path, checkpointFile)
}

#' CoGAPS with file input for matrix
#' @export
#'
#' @param D file path for data matrix
#' @return list with A and P matrix estimates
#' @examples
#'  file <- system.file("extdata/GIST.mtx", package="CoGAPS")
#'  CoGAPSFromFile(file)
CoGAPSFromFile <- function(D)
{
    cogapsFromFile_cpp(D)
}

#' Display Information About Package Compilation
#'
#' @details displays information about how the package was compiled, i.e. which
#'  compiler/version was used, which compile time options were enabled, etc...
#' @return display builds information
#' @examples
#'  CoGAPS::buildReport()
#' @export
buildReport <- function()
{
    getBuildReport_cpp()
}

#' Backwards Compatibility with v2
#'
#' @param D data matrix
#' @param S uncertainty matrix
#' @param ABins unused
#' @param PBins unused
#' @param simulation_id unused
#' @param nOutR number of output messages
#' @param output_atomic unused
#' @param fixedBinProbs unused
#' @param fixedDomain unused
#' @param sampleSnapshots indicates if snapshots should be made
#' @param numSnapshots how many snapshots to take
#' @param nMaxA unused
#' @param nMaxP unused
#' @param max_gibbmass_paraA limit truncated normal to max size
#' @param max_gibbmass_paraP limit truncated normal to max size
#' @return list with A and P matrix estimates
#' @importFrom methods new
#' @inheritParams CoGAPS
#' @examples
#' data(SimpSim)
#' result <- gapsRun(SimpSim.D, SimpSim.S, nFactor=3)
#' @export
gapsRun <- function(D, S, ABins=data.frame(), PBins=data.frame(), nFactor=7,
simulation_id="simulation", nEquil=1000, nSample=1000, nOutR=1000,
output_atomic=FALSE, fixedBinProbs=FALSE, fixedDomain="N", sampleSnapshots=TRUE,
numSnapshots=100, alphaA=0.01, nMaxA=100000, max_gibbmass_paraA=100.0,
alphaP=0.01, nMaxP=100000, max_gibbmass_paraP=100.0, seed=-1, messages=TRUE)
{
    warning('gapsRun is deprecated with v3.0, use CoGAPS')
    CoGAPS(D, S, nFactor=nFactor, nEquil=nEquil, nSample=nSample,
        nOutputs=nOutR, nSnapshots=ifelse(sampleSnapshots,numSnapshots,0),
        alphaA=alphaA, alphaP=alphaP, maxGibbmassA=max_gibbmass_paraA,
        messages=messages, maxGibbmassP=max_gibbmass_paraP, seed=seed)
}

#' Backwards Compatibility with v2
#'
#' @param D data matrix
#' @param S uncertainty matrix
#' @param FP data.frame with rows giving fixed patterns for P
#' @param fixedMatrix unused
#' @param ... v2 style parameters
#' @return list with A and P matrix estimates
#' @importFrom methods new
#' @inheritParams gapsRun
#' @examples
#' data(SimpSim)
#' nC <- ncol(SimpSim.D)
#' patterns <- matrix(1:nC/nC, nrow=1, ncol=nC)
#' result <- gapsMapRun(SimpSim.D, SimpSim.S, FP=patterns, nFactor=3)
#' @export
gapsMapRun <- function(D, S, FP, ABins=data.frame(), PBins=data.frame(),
nFactor=5, simulation_id="simulation", nEquil=1000, nSample=1000, nOutR=1000,
output_atomic=FALSE, fixedMatrix="P", fixedBinProbs=FALSE, fixedDomain="N",
sampleSnapshots=TRUE, numSnapshots=100, alphaA=0.01, nMaxA=100000,
max_gibbmass_paraA=100.0, alphaP=0.01, nMaxP=100000, max_gibbmass_paraP=100.0,
seed=-1, messages=TRUE)
{
    warning('gapsMapRun is deprecated with v3.0, use CoGaps')
    CoGAPS(D, S, nFactor=nFactor, nEquil=nEquil, nSample=nSample,
        nOutputs=nOutR, nSnapshots=ifelse(sampleSnapshots,numSnapshots,0),
        alphaA=alphaA, alphaP=alphaP, maxGibbmassA=max_gibbmass_paraA,
        messages=messages, maxGibbmassP=max_gibbmass_paraP, seed=seed,
        whichMatrixFixed='P', fixedPatterns=as.matrix(FP))
}

# helper function for backwards compatibility
v2CoGAPS <- function(result, ...)
{
    if (!is.null(list(...)$GStoGenes))
    {
        warning('GStoGenes is deprecated with v3.0, see CoGAPS documentation')
        if (is.null(list(...)$plot) | list(...)$plot)
        {
            plotGAPS(result$Amean, result$Pmean)
        }
        if (is.null(list(...)$nPerm))
        {
            nPerm <- 500
        }
        else
        {
            nPerm <- list(...)$nPerm
        }
        GSP <- calcCoGAPSStat(result$Amean, result$Asd, list(...)$GStoGenes,
            nPerm)
        result <- list(meanChi2=result$meanChi2, Amean=result$Amean,
            Asd=result$Asd, Pmean=result$Pmean, Psd=result$Psd,
            GSUpreg=GSP$GSUpreg, GSDownreg=GSP$GSDownreg, GSActEst=GSP$GSActEst)
    }
    return(result)
}
