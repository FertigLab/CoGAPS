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
#' @param ... additional parameters to be fed into \code{gapsRun} and \code{gapsMapRun}
#' @return list of A and P estimates
#' @seealso \code{\link{gapsRun}}, \code{\link{patternMatch4Parallel}}, and \code{\link{gapsMapRun}}
#' @examples
#' data(SimpSim)
#' sim_name <- "example"
#' createGWCoGAPSSets(SimpSim.D, SimpSim.S, nSets=2, sim_name)
#' result <- GWCoGAPS(sim_name, nFactor=3, nEquil=1000, nSample=1000)
#' @export
GWCoGAPS <- function(simulationName, nFactor, nCores=NA, cut=NA, minNS=NA, ...)
{
    if (!is.null(list(...)$checkpointFile))
        stop("checkpoint file name automatically set in GWCoGAPS - don't pass this parameter")
    allDataSets <- preInitialPhase(simulationName, nCores)
    initialResult <- runInitialPhase(simulationName, allDataSets, nFactor, ...)
    matchedPatternSets <- postInitialPhase(initialResult, length(allDataSets), cut, minNS)
    save(matchedPatternSets, file=paste(simulationName, "_matched_ps.RData", sep=""))
    finalResult <- runFinalPhase(simulationName, allDataSets, matchedPatternSets, ...)
    return(postFinalPhase(finalResult, matchedPatternSets))
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
        finalResult <- runFinalPhase(simulationName, allDataSets, matchedPatternSets, ...)
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
        finalResult <- runFinalPhase(simulationName, allDataSets, matchedPatternSets, ...)
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
        nCores <- length(allDataSets)
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
        cut <- nFactor
    return(patternMatch4Parallel(Ptot=BySet$P, nP=nFactor, nSets=nSets, cnt=cut,
        minNS=minNS, bySet=TRUE))
}

runFinalPhase <- function(simulationName, allDataSets, matchedPatternSets, ...)
{
    # generate seeds for parallelization
    nut <- generateSeeds(chains=length(allDataSets), seed=-1)

    # final number of factors
    nFactorFinal <- nrow(matchedPatternSets[[1]])

    # run fixed CoGAPS
    finalResult <- foreach(i=1:length(allDataSets)) %dopar%
    {
        # load data set and shift values so gene minimum is zero
        load(allDataSets[[i]])
        sampleD <- sweep(sampleD, 1, apply(sampleD, 1, function(x) pmin(0,min(x))))

        # run CoGAPS with fixed patterns
        cptFileName <- paste(simulationName, "_final_cpt_", i, ".out", sep="")
        CoGAPS(sampleD, sampleS, fixedPatterns=matchedPatternSets[[1]],
            nFactor=nFactorFinal, seed=nut[i], checkpointFile=cptFileName,
            whichMatrixFixed='P', ...)

    }
    return(finalResult)
}

postFinalPhase <- function(finalResult, matchedPatternSets)
{
    Aresult <- postFixed4Parallel(finalResult, matchedPatternSets[[1]])
    finalResult <- list("Amean"=Aresult$A, "Asd"=Aresult$Asd,
        "Pmean"=matchedPatternSets, "PbySet"=matchedPatternSets[["PBySet"]])
    class(finalResult) <- append(class(finalResult), "CoGAPS")
    return(finalResult)
}

