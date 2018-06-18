#' CoGAPS Matrix Factorization Algorithm
#' @export
#' @docType methods
#' @rdname CoGAPS-methods
#' 
#' @details calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix
#' @param D data
#' @param S uncertainty
#' @param params object of type CoGapsParams
#' @param ... keeps backwards compatibility with arguments from older versions
#' @return object of type CoGapsResult
#' @importFrom methods new
#' @examples
#' data(SimpSim)
<<<<<<< HEAD
#' result <- CoGAPS(SimpSim.D, SimpSim.S, nFactor=3, nOutputs=250)
#' @export
CoGAPS <- function(D, S, CoGAPSParams, GapsReturn ...)
{
    #returns a default uncertainty matrix of 0.1*D dataset if it's greater than 0.1
    if(is.null(S))
    S <- pmax(0.1*D, 0.1)
    
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
=======
#' result <- CoGAPS(SimpSim.D, nFactor=3, nOutputs=250)
setGeneric('CoGAPS', function(D, S=NULL, params, ...)
    {standardGeneric('RVsharing')})
>>>>>>> d4aca53ab9ee4ab023ea63049bda5a25e57b344b

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(D='matrix', params='CoGapsParams'),
function(D, S=NULL, params, ...)
{
    # default to 10% of signal if no uncertainty matrix passed in by user
    if (is.null(S))
    {
        S <- pmax(0.1 * D, 0.1)
    }

    return(RunCoGAPS(D, S, params))
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(D='data.frame', params='CoGapsParams'),
function(D, S=NULL, params, ...)
{
   
    D = data.matrix(D)
    
    if (is.null(S))
      S <- pmax(0.1 * D, 0.1)
    
    return(RunCoGAPS(D, S, params))
    # TODO(Hyejune) convert data.frame to matrix and set default value for
    # S matrix if it's null - then call RunCoGAPS
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(D='SummarizedExperiment', params='CoGapsParams'),
function(D, S=NULL, params, ...)
{
  
    D = assay(D, "counts")
    if (is.null(S))
      S <- pmax(0.1 * D, 0.1)
    return(RunCoGAPS(D, S, params))
    # TODO(Hyejune) extract count matrix from SE object and set default value
    # S matrix if it's null - then call RunCoGAPS
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(D='SingleCellExperiment', params='CoGapsParams'),
function(D, S=NULL, params, ...)
{
   
    D = assay(D, "counts")
    if (is.null(S))
      S <- pmax(0.1 * D, 0.1)
    return(RunCoGAPS(D, S, params))
    # TODO(Hyejune) extract count matrix from SCE object and set default value
    # S matrix if it's null - then call RunCoGAPS
})

RunCoGAPS <- function(D, S, params)
{
    thresholdEnum <- c("unique", "cut")
    result <- cogaps_cpp(D, S, params@nFactor, params@nIter, params@nIter / 10,
        params@nIter, params@outputFrequency, params@alphaA, params@alphaP,
        params@maxGibbmassA, params@maxGibbmassP, params@seed, params@messages,
        params@singleCellRNASeq, params@whichMatrixFixed, params@fixedMatrix,
        params@checkpointInterval, params@checkpointFile, 0, 0, params@nCores)

    return(new('CoGapsResult', result$Amean, result$Asd, result$Pmean,
        result$Psd))
}

# Add merge function so the results from CoGAPS can be incorporated 
# into a SingleCellExperiment
MergeResults <- function(result)
{  
  mergedResults <- merge(result$Amean, return$Asd, return$Pmean, return$Psd)
  resultsSCE <- SingleCellExperiment(assays=list(normcounts=cbind(result$Amean, result$Asd, result$Pmean, result$Psd))
  show(resultsSCE)
}
  
  
# TODO(Tom) handle instance where function is called with deprecated parametres