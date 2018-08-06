#' CoGAPS Distributed Matrix Factorization Algorithm
#'
#' @description runs CoGAPS over subsets of the data and stitches the results
#' back together
#' @details For file types CoGAPS supports csv, tsv, and mtx
#' @param data File name or R object (see details for supported types)
#' @param allParams list of all parameters used in computation
#' @param uncertainty uncertainty matrix (same supported types as data)
#' @return list
#' @importFrom BiocParallel bplapply MulticoreParam
distributedCogaps <- function(data, allParams, uncertainty)
{
    FUN <- function(index, sets, data, allParams, uncertainty, fixedMatrix=NULL)
    {
        internal <- ifelse(is(data, "character"), cogaps_cpp_from_file, cogaps_cpp)
        raw <- internal(data, allParams, uncertainty, sets[[index]],
            fixedMatrix, index == 1)
        new("CogapsResult", Amean=raw$Amean, Asd=raw$Asd, Pmean=raw$Pmean,
            Psd=raw$Psd, seed=raw$seed, meanChiSq=raw$meanChiSq)
    }

    # randomly sample either rows or columns into subsets to break the data up
    set.seed(allParams$gaps@seed)
    sets <- createSets(data, allParams)
    if (is.null(allParams$bpBackend))
        allParams$bpBackend <- BiocParallel::MulticoreParam(workers=length(sets))
    
    # run Cogaps normally on each subset of the data
    if (allParams$messages)
        message("Running Across Subsets...")
    initialResult <- bplapply(1:length(sets), FUN, BPPARAM=allParams$bpBackend,
        sets=sets, data=data, allParams=allParams, uncertainty=uncertainty)

    # allow user to save intermediate results
    if (allParams$saveUnmatched)
    {
        if (allParams$gaps@distributed == "genome-wide")
            unmatchedPatterns <- lapply(initialResult, function(x) x@sampleFactors)
        else
            unmatchedPatterns <- lapply(initialResult, function(x) x@featureLoadings)
        filename <- paste("unmatched_patterns_", allParams$gaps@nPatterns, sep="")
        save(sets, unmatchedPatterns, file=paste(filename, "RData", sep="."))
    }

    # match patterns in either A or P matrix
    if (allParams$messages)
        message("Matching Patterns Across Subsets...")
    consensusMatrix <- findConsensusMatrix(initialResult, allParams)
    allParams$gaps@nPatterns <- ncol(consensusMatrix)

    # set fixed matrix
    allParams$whichMatrixFixed <- ifelse(allParams$gaps@distributed
        == "genome-wide", "P", "A")

    # run final phase with fixed matrix
    if (allParams$messages)
        message("Running Final Stage...")
    finalResult <- bplapply(1:length(sets), FUN, BPPARAM=allParams$bpBackend,
        sets=sets, data=data, allParams=allParams, uncertainty=uncertainty,
        fixedMatrix=consensusMatrix)

    # get result 
    return(stitchTogether(finalResult, allParams))
}

#' partition genes/samples into subsets
#' @description either genes or samples or partitioned depending on the type
#' of distributed CoGAPS (i.e. genome-wide or single-cell)
#' @param data either file name or matrix
#' @param allParams list of all CoGAPS parameters
#' @return list of sorted subsets of either genes or samples
createSets <- function(data, allParams)
{
    total <- ifelse(xor(allParams$transposeData, allParams$gaps@distributed == "genome-wide"),
        nrow_helper(data), ncol_helper(data))
    setSize <- floor(total / allParams$gaps@nSets)

    sets <- list()
    remaining <- 1:total
    for (n in 1:(allParams$gaps@nSets - 1))
    {
        selected <- sample(remaining, setSize, replace=FALSE)
        sets[[n]] <- sort(selected)
        remaining <- setdiff(remaining, selected)
    }
    sets[[allParams$gaps@nSets]] <- sort(remaining)
    return(sets)
}

#' find the consensus pattern matrix across all subsets
#' @param result list of CogapsResult object from all runs across subsets
#' @param allParams list of all CoGAPS parameters
#' @return matrix of consensus patterns
findConsensusMatrix <- function(result, allParams)
{
    if (allParams$gaps@distributed == "genome-wide")
        patterns <- do.call(cbind, lapply(result, function(x) x@sampleFactors))
    else
        patterns <- do.call(cbind, lapply(result, function(x) x@featureLoadings))

    comb <- expand.grid(1:allParams$gaps@nSets, 1:allParams$gaps@nPatterns)
    colnames(patterns) <- paste(comb[,1], comb[,2], sep=".")
    return(patternMatch(patterns, allParams))
}

#' Match Patterns Across Multiple Runs
#' @param allPatterns matrix of patterns stored in the columns
#' @param allParams list of all CoGAPS parameters
#' @return a matrix of consensus patterns
#' @importFrom stats weighted.mean
patternMatch <- function(allPatterns, allParams)
{
    PatsByClust <- corcut(allPatterns, allParams)

    # split by maxNS
    indx <- which(sapply(PatsByClust, function(x) ncol(x) > allParams$gaps@maxNS))
    while (length(indx) > 0)
    { 
        allParams$gaps@cut <- 2
        internalPatsByClust <- corcut(PatsByClust[[indx[1]]], allParams)

        PatsByClust[[indx[1]]] <- internalPatsByClust[[1]]
        if (length(internalPatsByClust) > 1)
        {
            PatsByClust <- append(PatsByClust, internalPatsByClust[2])
        }
        indx <- which(sapply(PatsByClust, function(x) ncol(x) > allParams$gaps@maxNS))
    }

    # create matrix of mean patterns - weighted by coefficient of determination
    PatsByCDSWavg <- sapply(PatsByClust, function(clust)
        apply(clust, 1, function(row) weighted.mean(row, correlationToMeanPattern(clust)^3)))
    colnames(PatsByCDSWavg) <- paste("Pattern", 1:length(PatsByClust))

    # scale
    return(apply(PatsByCDSWavg, 2, function(col) col / max(col)))
}

correlationToMeanPattern <- function(cluster)
{
    meanPat <- rowMeans(cluster)
    sapply(1:ncol(cluster), function(j) round(cor(x=cluster[,j], y=meanPat), 3))
}

#' cluster patterns together
#' @param allPatterns matrix of all patterns across subsets
#' @param allParams list of all CoGAPS parameters
#' @return patterns listed by which cluster they belong to
#' @importFrom cluster agnes
#' @importFrom stats cutree as.hclust cor
corcut <- function(allPatterns, allParams)
{
    corr.dist <- cor(allPatterns)
    corr.dist <- 1 - corr.dist

    clust <- cluster::agnes(x=corr.dist, diss=TRUE, "complete")
    patternIds <- stats::cutree(stats::as.hclust(clust), k=allParams$gaps@cut)

    PatsByClust <- list()
    for (cluster in unique(patternIds))
    {
        if (sum(patternIds==cluster) >= allParams$gaps@minNS)
        {
            clusterPats <- allPatterns[,patternIds==cluster]
            PatsByClust[[as.character(cluster)]] <- clusterPats
        }
    }
    return(PatsByClust)
}

#' concatenate final results across subsets
#' @param result list of CogapsResult object from all runs across subsets
#' @param allParams list of all CoGAPS parameters
#' @return list with all CoGAPS output
stitchTogether <- function(result, allParams)
{
    if (allParams$gaps@distributed == "genome-wide")
    {
        consensus <- result[[1]]@sampleFactors
        Amean <- do.call(rbind, lapply(result, function(x) x@featureLoadings))
        Asd   <- do.call(rbind, lapply(result, function(x) x@featureStdDev))
        Pmean <- consensus
        Psd   <- matrix(0, nrow=nrow(consensus), ncol=ncol(consensus))
    }
    else
    {
        consensus <- result[[1]]@featureLoadings
        Amean <- consensus
        Asd   <- matrix(0, nrow=nrow(consensus), ncol=ncol(consensus))
        Pmean <- do.call(rbind, lapply(result, function(x) x@sampleFactors))
        Psd   <- do.call(rbind, lapply(result, function(x) x@sampleStdDev))
    }

    return(list("Amean"=Amean, "Asd"=Asd, "Pmean"=Pmean, "Psd"=Psd,
        "seed"=allParams$gaps@seed,
        "meanChiSq"=sum(sapply(result, function(r) r@metadata$meanChiSq))))
}

