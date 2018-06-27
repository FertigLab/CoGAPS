#' CoGAPS Matrix Factorization Algorithm
#' @export
#' @docType methods
#' @rdname CoGAPS-methods
#' 
#' @details calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix
#' @param Data data
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
CoGAPS <- function(Data, S, CoGAPSParams, GapsReturn ...)
{
    #returns a default uncertainty matrix of 0.1*Data dataset if it's greater than 0.1
    if(is.null(S))
    S <- pmax(0.1*Data, 0.1)
    
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
        Data <- oldArgs$data
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
setGeneric('CoGAPS', function(Data, S=NULL, params, ...)
    {standardGeneric('RVsharing')})
>>>>>>> d4aca53ab9ee4ab023ea63049bda5a25e57b344b

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(Data='matrix', params='CoGapsParams'),
function(Data, S=NULL, params, ...)
{
    # default to 10% of signal if no uncertainty matrix passed in by user
    if (is.null(S))
    {
        S <- pmax(0.1 * Data, 0.1)
    }

    return(RunCoGAPS(Data, S, params))
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(Data='data.frame', params='CoGapsParams'),
function(Data, S=NULL, params, ...)
{
   
    countMatrix = data.matrix(Data)
    
    if (is.null(S))
      S <- pmax(0.1 * Data, 0.1)
    
    return(RunCoGAPS(countMatrix, S, params))
    # TODO(Hyejune) convert data.frame to matrix and set default value for
    # S matrix if it's null - then call RunCoGAPS
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(Data='SummarizedExperiment', params='CoGapsParams'),
function(Data, S=NULL, params, ...)
{
  
    countMatrix = assay(Data, "counts")
    if (is.null(S))
      S <- pmax(0.1 * Data, 0.1)
    return(RunCoGAPS(countMatrix, S, params))
    # TODO(Hyejune) extract count matrix from SE object and set default value
    # S matrix if it's null - then call RunCoGAPS
})

#' @rdname CoGAPS-methods
#' @aliases CoGAPS
setMethod('CoGAPS', signature(Data='SingleCellExperiment', params='CoGapsParams'),
function(Data, S=NULL, params, ...)
{
   
    countMatrix = assay(Data, "counts")
    if (is.null(S))
      S <- pmax(0.1 * Data, 0.1)
    result <- RunCoGAPS(countMatrix, S, params)
    
    Data@reducedDims <- SimpleList(Amean=result$Amean, Pmean=result$Pmean)
    return(Data, Data@reducedDims)
    
    # TODO(Hyejune) extract count matrix from SCE object and set default value
    # S matrix if it's null - then call RunCoGAPS
})

RunCoGAPS <- function(Data, S, params)
{
    thresholdEnum <- c("unique", "cut")
    result <- cogaps_cpp(Data, S, params@nFactor, params@nIter, params@nIter / 10,
        params@nIter, params@outputFrequency, params@alphaA, params@alphaP,
        params@maxGibbmassA, params@maxGibbmassP, params@seed, params@messages,
        params@singleCellRNASeq, params@whichMatrixFixed, params@fixedMatrix,
        params@checkpointInterval, params@checkpointFile, 0, 0, params@nCores)

    return(new('CoGapsResult', result$Amean, result$Asd, result$Pmean,
        result$Psd))
}

  
# TODO(Tom) handle instance where function is called with deprecated parametres

