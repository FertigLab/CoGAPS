#' CoGAPS Distributed Matrix Factorization Algorithm
#' @keywords internal
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
    FUN <- function(index, sets, data, allParams, uncertainty, geneNames,
    sampleNames, fixedMatrix=NULL)
    {
        if (allParams$gaps@distributed == "genome-wide")
            geneNames <- geneNames[sets[[index]]]
        else
            sampleNames <- sampleNames[sets[[index]]]

        internal <- ifelse(is(data, "character"), cogaps_cpp_from_file, cogaps_cpp)
        raw <- internal(data, allParams, uncertainty, sets[[index]],
            fixedMatrix, index == 1)
        new("CogapsResult", Amean=raw$Amean, Asd=raw$Asd, Pmean=raw$Pmean,
            Psd=raw$Psd, meanChiSq=raw$meanChiSq, geneNames=geneNames,
            sampleNames=sampleNames)
    }

    # randomly sample either rows or columns into subsets to break the data up
    set.seed(allParams$gaps@seed)
    sets <- createSets(data, allParams)
    if (is.null(allParams$bpBackend))
        allParams$bpBackend <- BiocParallel::MulticoreParam(workers=length(sets))
    
    if (is.null(allParams$matchedPatterns))
    {
        # run Cogaps normally on each subset of the data
        if (allParams$messages)
            cat("Running Across Subsets...\n\n")
        initialResult <- bplapply(1:length(sets), FUN, BPPARAM=allParams$bpBackend,
            sets=sets, data=data, allParams=allParams, uncertainty=uncertainty,
            geneNames=allParams$geneNames, sampleNames=allParams$sampleNames)

        # get all unmatched patterns
        if (allParams$gaps@distributed == "genome-wide")
            unmatchedPatterns <- lapply(initialResult, function(x) x@sampleFactors)
        else
            unmatchedPatterns <- lapply(initialResult, function(x) x@featureLoadings)

        # match patterns in either A or P matrix
        if (allParams$messages)
            cat("\nMatching Patterns Across Subsets...\n")
        matchedPatterns <- findConsensusMatrix(unmatchedPatterns, allParams)
        allParams$gaps@nPatterns <- ncol(matchedPatterns$consensus)

        # set fixed matrix
        allParams$whichMatrixFixed <- ifelse(allParams$gaps@distributed
            == "genome-wide", "P", "A")
    }
    else
    {
        matchedPatterns <- list(consensus=allParams$matchedPatterns)
        allParams$gaps@nPatterns <- ncol(matchedPatterns$consensus)
        allParams$whichMatrixFixed <- ifelse(allParams$gaps@distributed
            == "genome-wide", "P", "A")
    }
        
    # run final phase with fixed matrix
    if (allParams$messages)
        cat("Running Final Stage...\n\n")
    finalResult <- bplapply(1:length(sets), FUN, BPPARAM=allParams$bpBackend,
        sets=sets, data=data, allParams=allParams, uncertainty=uncertainty,
        geneNames=allParams$geneNames, sampleNames=allParams$sampleNames,    
        fixedMatrix=matchedPatterns$consensus)

    # concatenate final result
    fullResult <- stitchTogether(finalResult, allParams)

    # add diagnostic information before returning
    if (is.null(allParams$matchedPatterns))
    {
        fullResult$diagnostics$unmatchedPatterns <- unmatchedPatterns
        fullResult$diagnostics$clusteredPatterns <- matchedPatterns$clusteredPatterns
        fullResult$diagnostics$CorrToMeanPattern <- lapply(matchedPatterns$clusteredPatterns, corrToMeanPattern)
    }

    if (allParams$gaps@distributed == "genome-wide")
        allNames <- allParams$geneNames
    else
        allNames <- allParams$sampleNames
    fullResult$diagnostics$subsets <- lapply(sets, function(s) allNames[s])

    # rename genes/samples if dimension was subsetted incompletely
    allUsedIndices <- sort(unlist(sets))
    if (allParams$gaps@distributed == "genome-wide")
    {
        fullResult$geneNames <- allParams$geneNames[allUsedIndices]
        fullResult$sampleNames <- allParams$sampleNames
    }
    else
    {
        fullResult$geneNames <- allParams$geneNames
        fullResult$sampleNames <- allParams$sampleNames[allUsedIndices]
    }
    return(fullResult)
}


#' find the consensus pattern matrix across all subsets
#' @keywords internal
#'
#' @param unmatchedPatterns list of all unmatched pattern matrices from initial
#' run of CoGAPS
#' @param allParams list of all CoGAPS parameters
#' @return matrix of consensus patterns
findConsensusMatrix <- function(unmatchedPatterns, allParams)
{
    allPatterns <- do.call(cbind, unmatchedPatterns)
    comb <- expand.grid(1:allParams$gaps@nSets, 1:allParams$gaps@nPatterns)
    colnames(allPatterns) <- paste(comb[,1], comb[,2], sep=".")
    return(patternMatch(allPatterns, allParams))
}

#' Match Patterns Across Multiple Runs
#' @keywords internal
#'
#' @param allPatterns matrix of patterns stored in the columns
#' @param allParams list of all CoGAPS parameters
#' @return a matrix of consensus patterns
#' @importFrom stats weighted.mean
patternMatch <- function(allPatterns, allParams)
{
    # cluster patterns
    clusters <- corcut(allPatterns, allParams$gaps@cut, allParams$gaps@minNS)

    # function to split a cluster in two (might fail to do so)
    splitCluster <- function(list, index, minNS)
    {
        split <- corcut(list[[index]], 2, minNS)
        list[[index]] <- split[[1]]
        if (length(split) > 1)
            list <- append(list, split[2])
        return(list)
    }

    # split large clusters into two
    tooLarge <- function(x) ncol(x) > allParams$gaps@maxNS
    indx <- which(sapply(clusters, tooLarge))
    while (length(indx) > 0)
    {
        clusters <- splitCluster(clusters, indx[1], allParams$gaps@minNS)
        indx <- which(sapply(clusters, tooLarge))
    }
    names(clusters) <- as.character(1:length(clusters))

    # create matrix of mean patterns - weighted by correlation to mean pattern
    meanPatterns <- sapply(clusters, function(clust) apply(clust, 1,
        function(row)  weighted.mean(row, corrToMeanPattern(clust)^3)))
    colnames(meanPatterns) <- paste("Pattern", 1:length(clusters))

    # returned patterns after scaling max to 1
    return(list("clusteredPatterns"=clusters,
        "consensus"=apply(meanPatterns, 2, function(col) col / max(col))))
}

#' calculate correlation of each pattern in a cluster to the cluster mean
#' @keywords internal
corrToMeanPattern <- function(cluster)
{
    meanPat <- rowMeans(cluster)
    sapply(1:ncol(cluster), function(j) round(cor(x=cluster[,j], y=meanPat), 3))
}

#' cluster patterns together
#' @keywords internal
#'
#' @param allPatterns matrix of all patterns across subsets
#' @param cut number of branches at which to cut dendrogram
#' @param minNS minimum of individual set contributions a cluster must contain
#' @return patterns listed by which cluster they belong to
#' @importFrom cluster agnes
#' @importFrom stats cutree as.hclust cor
corcut <- function(allPatterns, cut, minNS)
{
    corr.dist <- cor(allPatterns)
    corr.dist <- 1 - corr.dist

    clusterSummary <- cluster::agnes(x=corr.dist, diss=TRUE, "complete")
    clusterIds <- stats::cutree(stats::as.hclust(clusterSummary), k=cut)

    clusters <- list()
    for (id in unique(clusterIds))
    {
        if (sum(clusterIds==id) >= minNS)
            clusters[[as.character(id)]] <- allPatterns[,clusterIds==id,drop=FALSE]
    }
    return(clusters)
}

#' concatenate final results across subsets
#' @keywords internal
#'
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

