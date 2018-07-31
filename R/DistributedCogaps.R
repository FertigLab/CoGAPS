#' CoGAPS Distributed Matrix Factorization Algorithm
#'
#' @description runs CoGAPS over subsets of the data and stitches the results
#' back together
#' @details For file types CoGAPS supports csv, tsv, and mtx
#' @param data File name or R object (see details for supported types)
#' @param allParams list of all parameters used in computation
#' @param uncertainty uncertainty matrix (same supported types as data)
#' @return list
#' @importFrom BiocParallel bplapply SnowParam
distributedCogaps <- function(data, allParams, uncertainty)
{
    FUN <- function(set, data, allParams, uncertainty, fixedMatrix=NULL)
    {
        internal <- ifelse(is(data, "character"), cogaps_cpp_from_file, cogaps_cpp)
        raw <- internal(data, allParams, uncertainty, set, fixedMatrix)
        new("CogapsResult", Amean=raw$Amean, Asd=raw$Asd, Pmean=raw$Pmean,
            Psd=raw$Psd, seed=raw$seed, meanChiSq=raw$meanChiSq)    
    }

    # randomly sample either rows or columns into subsets to break the data up
    set.seed(allParams$gaps@seed)
    sets <- createSets(data, allParams)
    snow <- SnowParam(workers=length(sets), type="SOCK")

    # run Cogaps normally on each subset of the data
    initialResult <- bplapply(sets, FUN, BPPARAM=snow, data=data,
        allParams=allParams, uncertainty=uncertainty)

    # match patterns in either A or P matrix
    consensusMatrix <- findConsensusMatrix(initialResult, allParams)
    allParams$gaps@nPatterns <- ncol(consensusMatrix)

    # set fixed matrix
    allParams$whichMatrixFixed <- ifelse(allParams$gaps@distributed
        == "genome-wide", "P", "A")

    # ru final phase with fixed matrix
    finalResult <- bplapply(sets, FUN, data=data, BPPARAM=snow,
        allParams=allParams, uncertainty=uncertainty, fixedMatrix=consensusMatrix)

    # get result 
    return(stitchTogether(finalResult, allParams))
}

#' get number of rows from supported file name or matrix
#' @param data either a file name or a matrix
#' @return number of rows
#' @importFrom data.table fread
#' @importFrom tools file_ext
nrow_helper <- function(data)
{
    if (is(data, "character"))
    {
        return(switch(tools::file_ext(data),
            "csv" = nrow(data.table::fread(data, select=1)),
            "tsv" = nrow(data.table::fread(data, select=1)),
            "mtx" = as.numeric(data.table::fread(data, nrows=1, fill=TRUE)[1,1])
        ))
    }
    return(nrow(data))
}

#' get number of columns from supported file name or matrix
#' @param data either a file name or a matrix
#' @return number of columns
#' @importFrom data.table fread
#' @importFrom tools file_ext
ncol_helper <- function(data)
{
    if (is(data, "character"))
    {
        return(switch(tools::file_ext(data),
            "csv" = ncol(data.table::fread(data, nrows=1)) - 1,
            "tsv" = ncol(data.table::fread(data, nrows=1)) - 1,
            "mtx" = as.numeric(data.table::fread(data, nrows=1, fill=TRUE)[1,2])
        ))
    }
    return(ncol(data))
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
#' @importFrom stats weighted.mean
patternMatch <- function(allPatterns, allParams)
{
    cc <- corcut(allPatterns, allParams)

    # split by maxNS
    indx <- which(sapply(cc$PatsByClust, function(x) ncol(x) > allParams$gaps@maxNS))
    while (length(indx) > 0)
    { 
        allParams$gaps@cut <- 2
        icc <- corcut(cc$PatsByClust[[indx[1]]], allParams)

        cc$PatsByClust[[indx[1]]] <- icc$PatsByClust[[1]]
        cc$RtoMeanPattern[[indx[1]]] <- icc$RtoMeanPattern[[1]]
        if (length(icc$PatsByClust) > 1)
        {
            cc$PatsByClust <- append(cc$PatsByClust, icc$PatsByClust[2])
            cc$RtoMeanPattern <- append(cc$RtoMeanPattern, icc$RtoMeanPattern[2])
        }
        indx <- which(sapply(cc$PatsByClust, function(x) ncol(x) > allParams$gaps@maxNS))
    }

    # create matrix of mean patterns - weighted by coefficient of determination
    PatsByCDSWavg <- sapply(1:length(cc$PatsByClust), function(z)
        apply(cc$PatsByClust[[z]], 1, function(x) weighted.mean(x, (cc$RtoMeanPattern[[z]])^3)))
    colnames(PatsByCDSWavg) <- lapply(1:length(cc$PatsByClust), function(x) paste("Pattern", x))

    # scale
    return(apply(PatsByCDSWavg, 2, function(col) col / max(col)))
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

