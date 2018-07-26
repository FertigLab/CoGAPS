#' CoGAPS Distributed Matrix Factorization Algorithm
#'
#' @description runs CoGAPS over subsets of the data and stitches the results
#' back together
#' @details For file types CoGAPS supports csv, tsv, and mtx
#' @param data File name or R object (see details for supported types)
#' @param allParams list of all parameters used in computation
#' @param uncertainty uncertainty matrix (same supported types as data)
#' @return list
distributedCogaps <- function(data, allParams, uncertainty)
{
    # randomly sample either rows or columns into subsets to break the data up
    sets <- createSets(data, allParams)

    # run Cogaps normally on each subset of the data
    initialResult <- foreach(i=1:allParams$gaps@nSets) %dopar%
    {
        cogaps_cpp(data, allParams, uncertainty, sets[[i]])
    }

    # match patterns in either A or P matrix
    consensusMatrix <- findConsensusMatrix(initialResult, allParams)

    # set fixed matrix
    allParams$whichMatrixFixed <- ifelse(allParams$gaps@distributed
        == "genome-wide", "P", "A")

    # run all subsets with the same fixed matrix
    finalResult <- foreach(i=1:nSets) %dopar%
    {
        cogaps_cpp(data, allParams, uncertainty, sets[[i]], consensusMatrix)
    }

    # get result 
    resultList <- stitchTogether(finalResult, allParams)
    resultList$seed <- allParams$modelParams@seed
    resultList$meanChiSq <- sum(sapply(finalResult, function(r) r$meanChiSq))
    resultList$diagnostics <- list()
    return(resultList)
}

nrow_helper <- function(data)
{
    if (class(data) == "character")
    {
        switch(file_ext(data),
            "csv" = nrow(data.table::fread(data, select=1)),
            "tsv" = nrow(data.table::fread(data, select=1)),
            "mtx" = as.numeric(data.table::fread(data, nrows=1, fill=TRUE)[1,1])
        )
    }
    return(nrow(data))
}

ncol_helper <- function(data)
{
    if (class(data) == "character")
    {
        switch(file_ext(data),
            "csv" = ncol(data.table::fread(data, nrows=1)) - 1,
            "tsv" = ncol(data.table::fread(data, nrows=1)) - 1,
            "mtx" = as.numeric(data.table::fread(data, nrows=1, fill=TRUE)[1,2])
        )
    }
    return(ncol(data))
}        

createSets <- function(data, allParams)
{
    total <- ifelse(allParams$gaps@distributed == "genome-wide",
        nrow_helper(data), ncol_helper(data))
    setSize <- floor(total / allParams$gaps@nSets)

    sets <- list()
    remaining <- 1:total
    for (n in 1:(allParams$gaps@nSets - 1))
    {
        selected <- sample(remaining, setSize, replace=FALSE)
        sets[[n]] <- selected
        remaining <- setdiff(remaining, selected)
    }
    sets[[allParams$gaps@nSets]] <- remaining
    return(sets)
}

findConsensusMatrix <- function(result, allParams)
{
    if (allParams$gaps@distributed == "genome-wide")
        patterns <- do.call(cbind, lapply(result, function(x) x@sampleFactors))
    else
        patterns <- do.call(cbind, lapply(result, function(x) x@featureLoadings))

    comb <- expand.grid(1:allParams$gaps@nSets, 1:allParams$gaps@nPatterns)
    rownames(patterns) <- paste(comb[,1], comb[,2], sep=".")
    return(patternMatch(patterns, allParams))
}

#' Match Patterns Across Multiple Runs
#' @export
#'
#' @param Atot a matrix containing the total by set estimates of Pmean output from \code{reOrderBySet}
#' @param nSets number of parallel sets used to generate \code{Atot}
#' @param cnt  number of branches at which to cut dendrogram
#' @param minNS minimum of individual set contributions a cluster must contain
#' @param maxNS maximum of individual set contributions a cluster must contain
#' @param ignore.NA logical indicating whether or not to ignore NAs from potential over dimensionalization. Default is FALSE.
#' @param bySet logical indicating whether to return list of matched set solutions from \code{Atot}
#' @param ... additional parameters for \code{agnes}
#' @return a matrix of consensus patterns by samples. If \code{bySet=TRUE} then a list of the set contributions to each
#' consensus pattern is also returned.
patternMatch <- function(allPatterns, allParams)
{
    cc <- corcut(allPatterns, allParams)

    ### split by maxNS
    indx <- which(sapply(cc$PatsByClust, function(x) ncol(x) > allParams$gaps@maxNS))
    while (length(indx) > 0)
    { 
        allParams$gaps@cut <- 2
        icc <- corcut(cc$PatsByClust[[indx[1]]], allParams)

        cc$PatsByClust[[indx[1]]] <- icc$PatsByClust[[1]]
        cc$RtoMeanPattern[[indx[1]]] <- icc$RtoMeanPattern[[1]]
        cc$PatsByClust <- append(cc$PatsByClust, icc$PatsByClust[2])
        cc$RtoMeanPattern <- append(cc$RtoMeanPattern, icc$RtoMeanPattern[2])

        indx <- which(sapply(cc$PatsByClust, function(x) ncol(x) > allParams$gaps@maxNS))
    }

    # created matrix of mean patterns - weighted by R closeness of fit
    PatsByCDSWavg <- t(sapply(1:length(cc$PatsByClust), function(z)
        apply(cc$PatsByClust[[z]], 1, function(x) weighted.mean(x, (cc$RtoMeanPattern[[z]])^3))))
    rownames(PatsByCDSWavg) <- lapply(1:length(cc$PatsByClust), function(x) paste("Pattern", x))
    return(PatsByCDSWavg)

    # scale
    return(apply(PatsByCDSWavg, 1, function(row) row / max(row)))
}

corcut <- function(allPatterns, allParams)
{
    corr.dist <- cor(allPatterns)
    corr.dist <- 1 - corr.dist

    clust <- agnes(x=corr.dist, diss=TRUE, "complete")
    patternIds <- cutree(as.hclust(clust), k=allParams$gaps@cut)

    RtoMeanPattern <- list()
    PatsByClust <- list()
    for (cluster in unique(patternIds))
    {
        if (sum(patternIds==cluster) >= allParams$gaps@minNS)
        {
            clusterPats <- allPatterns[,patternIds==cluster]
            meanPat <- rowMeans(clusterPats)
            RtoMeanPattern[[as.character(cluster)]] <- sapply(1:ncol(clusterPats),
                function(j) round(cor(x=clusterPats[,j], y=meanPat), 3))
            PatsByClust[[as.character(cluster)]] <- clusterPats
        }
    }
    return(list("RtoMeanPattern"=RtoMeanPattern, "PatsByClust"=PatsByClust))
}

stitchTogether <- function(result, allParams)
{
    if (allParams$modelParams@distributed == "genome-wide")
    {
        consensus <- result[[1]]@Pmean
        return(list(
            "Amean" = do.call(rbind, lapply(result, function(x) x@featureLoadings)),
            "Asd"   = do.call(rbind, lapply(result, function(x) x@featureStdDev)),
            "Pmean" = consensus,
            "Psd"   = matrix(0, nrow=nrow(consensus), ncol=ncol(consensus))
        ))
    }
    else
    {
        consensus <- result[[1]]@Amean
        return(list(
            "Amean" = consensus,
            "Asd"   = matrix(0, nrow=nrow(consensus), ncol=ncol(consensus)),
            "Pmean" = do.call(rbind, lapply(result, function(x) x@sampleFactors)),
            "Psd"   = do.call(rbind, lapply(result, function(x) x@sampleStdDev))
        ))
    }
}

