CoGAPS <- function(D, S=NULL, nFactor=7, nEquil=1000, nSample=1000, nOutputs=1000,
                   nSnapshots=0, alphaA=0.01, alphaP=0.01, maxGibbmassA=100, maxGibbmassP=100,
                   seed=-1, messages=TRUE, singleCellRNASeq=FALSE, whichMatrixFixed='N',
                   fixedPatterns=matrix(0), checkpointInterval=0, 
                   checkpointFile="gaps_checkpoint.out", nCores=1, ...)
{
  #returns a default uncertainty matrix of 0.1*D dataset if it's greater than 0.1
  if(!is.null(S))
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
  if (seed < 0)
  {
    # TODO get time in milliseconds
    seed <- 0
  }
  
  # run algorithm with call to C++ code
  result <- cogaps_cpp(D, S, nFactor, nEquil, nEquil/10, nSample, nOutputs,
                       nSnapshots, alphaA, alphaP, maxGibbmassA, maxGibbmassP, seed, messages,
                       singleCellRNASeq, whichMatrixFixed, fixedPatterns, checkpointInterval,
                       checkpointFile, which(thresholdEnum==pumpThreshold), nPumpSamples,
                       nCores)
  
  # label matrices and return list
  patternNames <- paste('Patt', 1:nFactor, sep='')
  rownames(result$Amean) <- rownames(result$Asd) <- rownames(D)
  colnames(result$Amean) <- colnames(result$Asd) <- patternNames
  rownames(result$Pmean) <- rownames(result$Psd) <- patternNames
  colnames(result$Pmean) <- colnames(result$Psd) <- colnames(D)
  return(v2CoGAPS(result, ...)) # backwards compatible with v2
}
  
