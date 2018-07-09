#' scCoGAPS
#' @export
#'
#' @details calls the C++ MCMC code and performs Bayesian
#' matrix factorization returning the two matrices that reconstruct
#' the data matrix for whole genome data;
#' @param simulationName name of this simulation
#' @param nFactor number of patterns (basis vectors, metagenes), which must be
#' greater than or equal to the number of rows of FP
#' @param nCores number of cores for parallelization. If left to the default NA, nCores = nSets.
#' @param cut number of branches at which to cut dendrogram used in patternMatch4singleCell
#' @param minNS minimum of individual set contributions a cluster must contain
#' @param manualMatch logical indicating whether or not to stop after initial phase for manual pattern matching
#' @param consensusAs fixed pattern matrix to be used to ensure reciprocity of A weights accross sets 
#' @param ... additional parameters to be fed into \code{gapsRun} and \code{gapsMapRun}
#' @return list of A and P estimates
scCoGAPS <- function(simulationName, nFactor, nCores=NA, cut=NA, minNS=NA, manualMatch=FALSE, consensusAs=NULL, ...)
{
    if (!is.null(list(...)$checkpointFile))
    {
        stop("checkpoint file name automatically set in GWCoGAPS - don't pass this parameter")
    }

    if (is.null(consensusAs))
    {
        allDataSets <- sc_preInitialPhase(simulationName, nCores)
        initialResult <- sc_runInitialPhase(simulationName, allDataSets, nFactor, ...)
        if (manualMatch)
        {
            saveRDS(initialResult,file=paste(simulationName,"_initial.rds", sep=""))
            stop("Please provide concensus gene weights upon restarting.")
        }
        matchedAmplitudes <- sc_postInitialPhase(initialResult, length(allDataSets), cut, minNS)
        consensusAs <- matchedAmplitudes[[1]]
        save(consensusAs, file=paste(simulationName, "_matched_As.RData", sep=""))
    } 
    finalResult <- sc_runFinalPhase(simulationName, allDataSets, consensusAs, nCores, ...)
    return(sc_postFinalPhase(finalResult, consensusAs))
}

#' Restart a scCoGAPS run from a Checkpoint
#'
#' @inheritParams GWCoGAPS
#' @return list of A and P estimates
#' @importFrom utils file_test
scCoGapsFromCheckpoint <- function(simulationName, nCores, cut=NA, minNS=NA, ...)
{
    # find data files
    allDataSets <- sc_preInitialPhase(simulationName, nCores)

    # figure out phase from file signature
    initialCpts <- list.files(full.names=TRUE, pattern=paste(simulationName,
            "_initial_cpt_[0-9]+.out", sep=""))
    finalCpts <- list.files(full.names=TRUE, pattern=paste(simulationName,
            "_final_cpt_[0-9]+.out", sep=""))

    if (length(finalCpts))
    {
        finalResult <- foreach(i=1:length(allDataSets)) %dopar%
        {
            # load data set and shift values so gene minimum is zero
            load(allDataSets[[i]])
            sampleD <- sweep(sampleD, 1, apply(sampleD, 1, function(x)
                pmin(0,min(x))))
    
            # run CoGAPS with fixed patterns
            cptFileName <- paste(simulationName, "_final_cpt_", i, ".out",
                sep="")
            CoGapsFromCheckpoint(sampleD, sampleS, cptFileName)
        }
        load(paste(simulationName, "_matched_As.RData", sep=""))
    }
    else if (file_test("-f", paste(simulationName, "_matched_As.RData", sep="")))
    {
        load(paste(simulationName, "_matched_As.RData", sep=""))
        finalResult <- sc_runFinalPhase(simulationName, allDataSets,
            consensusAs, nCores, ...)
    }
    else if (length(initialCpts))
    {
        # initial phase - always needs to be run to get matchedAmplitudes
        initialResult <- foreach(i=1:length(allDataSets)) %dopar%
        {
            # load data set and shift values so gene minimum is zero
            load(allDataSets[[i]])
            sampleD <- sweep(sampleD, 1, apply(sampleD, 1, function(x)
                pmin(0,min(x))))

            # run CoGAPS from checkpoint
            cptFileName <- paste(simulationName, "_initial_cpt_", i, ".out",
                sep="")
            CoGapsFromCheckpoint(sampleD, sampleS, cptFileName)
        }
        matchedAmplitudes <- sc_postInitialPhase(initialResult,
            length(allDataSets), cut, minNS)
        consensusAs <- matchedAmplitudes[[1]]
        save(consensusAs, file=paste(simulationName, "_matched_As.RData",
            sep=""))
        finalResult <- sc_runFinalPhase(simulationName, allDataSets,
            consensusAs, nCores, ...)
    }
    else
    {
        stop("no checkpoint files found")
    }
    return(sc_postFinalPhase(finalResult, consensusAs))
}

sc_preInitialPhase <- function(simulationName, nCores)
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

sc_runInitialPhase <- function(simulationName, allDataSets, nFactor, ...)
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
            checkpointFile=cptFileName, singleCellRNASeq=TRUE, ...)
    }
    return(initialResult)
}

sc_postInitialPhase <- function(initialResult, nSets, cut, minNS)
{
    nFactor <- ncol(initialResult[[1]]$Amean)
    BySet <- reOrderBySet(AP=initialResult, nFactor=nFactor, nSets=nSets,match="A")

    #run postpattern match function
    if (is.na(cut))
    {
        cut <- nFactor
    }

    return(cellMatchR(Atot=BySet$A, nSets=nSets, cnt=cut, minNS=minNS,
        bySet=TRUE))
}

sc_runFinalPhase <- function(simulationName, allDataSets, consensusAs, nCores, ...)
{    
    if (length(dim(consensusAs)) != 2)
    {
        stop("consensusAs must be a matrix")
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
    nFactorFinal <- ncol(consensusAs)

    # run fixed CoGAPS
    finalResult <- foreach(i=1:length(allDataSets)) %dopar%
    {
        # load data set and shift values so gene minimum is zero
        load(allDataSets[[i]])
        sampleD <- sweep(sampleD, 1, apply(sampleD, 1, function(x) pmin(0,min(x))))

        # run CoGAPS with fixed patterns
        cptFileName <- paste(simulationName, "_final_cpt_", i, ".out", sep="")
        CoGAPS(sampleD, uncertainty=sampleS, fixedMatrix=consensusAs,
            nFactor=nFactorFinal, seed=nut[i], checkpointFile=cptFileName,
            whichMatrixFixed='A', singleCellRNASeq=TRUE, ...)

    }
    return(finalResult)
}

sc_postFinalPhase <- function(finalResult, consensusAs)
{
    Aresult <- postFixed4Parallel(finalResult, consensusAs, setMatrix="A")
    finalResult <- list("Pmean"=Aresult$P, "Psd"=Aresult$Psd,"Amean"=consensusAs)
    class(finalResult) <- append(class(finalResult), "CoGAPS")
    return(finalResult)
}
