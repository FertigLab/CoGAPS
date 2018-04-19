#' GWCoGAPS
#'
#' @details calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix for whole genome data;
#' @param simulationName name of this simulation
#' @param nFactor number of patterns (basis vectors, metagenes), which must be
#' greater than or equal to the number of rows of FP
#' @param nCores number of cores for parallelization. If left to the default NA, nCores = nSets.
#' @param cut number of branches at which to cut dendrogram used in patternMatch4Parallel
#' @param minNS minimum of individual set contributions a cluster must contain
#' @param manualMatch logical indicating whether or not to stop after initial phase for manual pattern matching
#' @param consensusPatterns fixed pattern matrix to be used to ensure reciprocity of A weights accross sets 
#' @param ... additional parameters to be fed into \code{gapsRun} and \code{gapsMapRun}
#' @return list of A and P estimates
#' @seealso \code{\link{gapsRun}}, \code{\link{patternMatch4Parallel}}, and \code{\link{gapsMapRun}}
#' @examples
#' data(SimpSim)
#' sim_name <- "example"
#' createGWCoGAPSSets(SimpSim.D, SimpSim.S, nSets=2, sim_name)
#' result <- GWCoGAPS(sim_name, nFactor=3, nEquil=200, nSample=200)
#' @export
GWCoGAPS <- function(simulationName, nFactor, nCores=NA, cut=NA, minNS=NA, manualMatch=FALSE, consensusPatterns=NULL, ...)
{
    if (!is.null(list(...)$checkpointFile))
    {
        stop("checkpoint file name automatically set in GWCoGAPS - don't pass this parameter")
    }

    if (is.null(consensusPatterns))
    {
        allDataSets <- preInitialPhase(simulationName, nCores)
        initialResult <- runInitialPhase(simulationName, allDataSets, nFactor, ...)
        if (manualMatch)
        {
            saveRDS(initialResult,file=paste(simulationName,"_initial.rds", sep=""))
            stop("Please provide consensus patterns upon restarting.")
        }
        matchedPatternSets <- postInitialPhase(initialResult, length(allDataSets), cut, minNS)
        save(matchedPatternSets, file=paste(simulationName, "_matched_ps.RData", sep=""))
        consensusPatterns<-matchedPatternSets[[1]]
    } 
    finalResult <- runFinalPhase(simulationName, allDataSets, consensusPatterns, nCores, ...)
    return(postFinalPhase(finalResult, consensusPatterns))
}

#' Restart a GWCoGaps Run from Checkpoint
#'
#' @inheritParams GWCoGAPS
#' @return list of A and P estimates
#' @importFrom utils file_test
#' @export
GWCoGapsFromCheckpoint <- function(simulationName, nCores=NA, cut=NA, minNS=NA, ...)
{
    # find data files
    allDataSets <- preInitialPhase(simulationName, nCores)

    # figure out phase from file signature
    initialCpts <- list.files(full.names=TRUE, pattern=paste(simulationName, "_initial_cpt_[0-9]+.out", sep=""))
    finalCpts <- list.files(full.names=TRUE, pattern=paste(simulationName, "_final_cpt_[0-9]+.out", sep=""))

    if (length(finalCpts))
    {
        finalResult <- foreach(i=1:length(allDataSets)) %dopar%
        {
            # load data set and shift values so gene minimum is zero
            load(allDataSets[[i]])
            sampleD <- sweep(sampleD, 1, apply(sampleD, 1, function(x) pmin(0,min(x))))
    
            # run CoGAPS with fixed patterns
            cptFileName <- paste(simulationName, "_final_cpt_", i, ".out", sep="")
            CoGapsFromCheckpoint(sampleD, sampleS, cptFileName)
        }
        load(paste(simulationName, "_matched_ps.RData", sep=""))
    }
    else if (file_test("-f", paste(simulationName, "_matched_ps.RData", sep="")))
    {
        load(paste(simulationName, "_matched_ps.RData", sep=""))
        consensusPatterns<-matchedPatternSets[[1]]
        finalResult <- runFinalPhase(simulationName, allDataSets, consensusPatterns, ...)
    }
    else if (length(initialCpts))
    {
        # initial phase - always needs to be run to get matchedPatternSets
        initialResult <- foreach(i=1:length(allDataSets)) %dopar%
        {
            # load data set and shift values so gene minimum is zero
            load(allDataSets[[i]])
            sampleD <- sweep(sampleD, 1, apply(sampleD, 1, function(x) pmin(0,min(x))))

            # run CoGAPS from checkpoint
            cptFileName <- paste(simulationName, "_initial_cpt_", i, ".out", sep="")
            CoGapsFromCheckpoint(sampleD, sampleS, cptFileName)
        }
        matchedPatternSets <- postInitialPhase(initialResult, length(allDataSets), cut, minNS)
        save(matchedPatternSets, file=paste(simulationName, "_matched_ps.RData", sep=""))
        consensusPatterns<-matchedPatternSets[[1]]
        finalResult <- runFinalPhase(simulationName, allDataSets, consensusPatterns, ...)
    }
    else
    {
        stop("no checkpoint files found")
    }
    return(postFinalPhase(finalResult, matchedPatternSets))
}

preInitialPhase <- function(simulationName, nCores)
{
    # find data files
    fileSig <- paste(simulationName, "_partition_[0-9]+.RData", sep="")
    allDataSets <- list.files(full.names=TRUE, pattern=fileSig)

    # establish the number of cores that we are able to use
    if (is.na(nCores))
    {
        nCores <- length(allDataSets)
    }
    registerDoParallel(cores=nCores)
    return(allDataSets)
}

runInitialPhase <- function(simulationName, allDataSets, nFactor, ...)
{
    #generate seeds for parallelization
    nut <- generateSeeds(chains=length(allDataSets), seed=-1)

    # run CoGAPS for each set
    initialResult <- foreach(i=1:length(allDataSets)) %dopar%
    {
        # load data set and shift values so gene minimum is zero
        load(allDataSets[[i]])
        sampleD <- sweep(sampleD, 1, apply(sampleD, 1, function(x) pmin(0,min(x))))

        # run CoGAPS without any fixed patterns
        cptFileName <- paste(simulationName, "_initial_cpt_", i, ".out", sep="")
        CoGAPS(sampleD, sampleS, nFactor=nFactor, seed=nut[i],
            checkpointFile=cptFileName, ...)
    }
    return(initialResult)
}

postInitialPhase <- function(initialResult, nSets, cut, minNS)
{
    nFactor <- ncol(initialResult[[1]]$Amean)
    BySet <- reOrderBySet(AP=initialResult, nFactor=nFactor, nSets=nSets)

    #run postpattern match function
    if (is.na(cut))
    {
        cut <- nFactor
    }

    return(patternMatch4Parallel(Ptot=BySet$P, nP=nFactor, nSets=nSets, cnt=cut,
        minNS=minNS, bySet=TRUE))
}

runFinalPhase <- function(simulationName, allDataSets, consensusPatterns, nCores, ...)
{    
    if (length(dim(consensusPatterns)) != 2)
    {
        stop("consensusPatterns must be a matrix")
    }

    # find data files if providing consensus patterns
    fileSig <- paste(simulationName, "_partition_[0-9]+.RData", sep="")
    allDataSets <- list.files(full.names=TRUE, pattern=fileSig)

    if (is.na(nCores))
    {
        nCores <- length(allDataSets)
    }
    registerDoParallel(cores=nCores)

    # generate seeds for parallelization
    nut <- generateSeeds(chains=length(allDataSets), seed=-1)

    # final number of factors
    nFactorFinal <- nrow(consensusPatterns)

    # run fixed CoGAPS
    finalResult <- foreach(i=1:length(allDataSets)) %dopar%
    {
        # load data set and shift values so gene minimum is zero
        load(allDataSets[[i]])
        sampleD <- sweep(sampleD, 1, apply(sampleD, 1, function(x) pmin(0,min(x))))

        # run CoGAPS with fixed patterns
        cptFileName <- paste(simulationName, "_final_cpt_", i, ".out", sep="")
        CoGAPS(sampleD, sampleS, fixedPatterns=consensusPatterns,
            nFactor=nFactorFinal, seed=nut[i], checkpointFile=cptFileName,
            whichMatrixFixed='P', ...)

    }
    return(finalResult)
}

postFinalPhase <- function(finalResult, consensusPatterns)
{
    Aresult <- postFixed4Parallel(finalResult, consensusPatterns)
    finalResult <- list("Amean"=Aresult$A, "Asd"=Aresult$Asd,"Pmean"=consensusPatterns)
    class(finalResult) <- append(class(finalResult), "CoGAPS")
    return(finalResult)
}

