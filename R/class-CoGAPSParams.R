---
  title: "class-CoGAPSParams"
author: "Hyejune Limb"
date: "6/11/2018"
output: html_document
---
  
  ```{r}

CoGAPSParams <- setClass("CoGAPSParams", slots = c(
  nFactor = "numeric",
  nEquil = "numeric",
  nSample = "numeric",
  nOutputs = "numeric",
  nSnapshots = "numeric",
  alphaA = "numeric",
  alphaP = "numeric",
  maxGibbmassA = "numeric",
  maxGibbMassP = "numeric",
  seed = "numeric",
  messages = "numeric",
  singleCellRNASeq = "logical",
  whichMatrixFixed = "character",
  fixedPatterns = "matrix",
  checkpointInterval = "numeric",
  checkpointFile = "character", 
  nCores = "numeric"
))


setMethod("initialize", "CoGAPSParams",
          function(.Object, ...)
          {
            .Object@nFactor <- 7
            .Object@nEquil <- 1000
            .Object@nSample <- 1000
            .Object@nOutputs <- 1000
            .Object@nSnapshots <- 0
            .Object@alphaA <- 0.01
            .Object@alphaP <- 0.01
            .Object@maxGibbmassA <- 100
            .Object@maxGibbmassP <- 100
            .Object@seed <- -1
            .Object@messages <- TRUE
            .Object@singleCellRNASeq <- FALSE
            .Object@whichMatrixFixed <- 'N'
            .Object@fixedPatterns <- matrix(0)
            .Object@checkpointInterval <- 0
            .Object@checkpointFile <- "gaps_checkpoint.out"
            .Object@nCores <- 1
            
            .Object <- callNextMethod(.Object, ...)
            .Object
          }
)

setValidity("CoGAPSParams",
            function(object)
            {
              if (object@nFactor < 0)
                stop('number of patterns must be non-negative')
              if (object@nEquil < 0)
                stop('number of iterations for burn-in must be non-negative')
              if (object@nSample < 0)
                stop('number of iterations for sampling must be non-negative')
              if (object@nOutputs < 0)
                stop('number of iterations to print status into R by iterations must be non-negative')
              if (object@nSnapshots < 0)
                stop('number of samples to capture must be non-negative')
              if (object@alphaA  < 0)
                stop('sparsity parameter for A domain must be non-negative')
              if (object@alphaP  < 0)
                stop('sparsity parameter for P domain must be non-negative')
              if (object@maxGibbmassA < 0)
                stop('limit must be non-negative')
              if (object@maxGibbmassP < 0)
                stop('limit must be non-negative')
              if (object@nCores < 0)
                stop('number of cpu cores to run in parallel order must be non-negative')
              
            }
)

```