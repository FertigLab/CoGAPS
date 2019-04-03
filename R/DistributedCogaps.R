#' make correct call to internal CoGAPS dispatch function, CoGAPS could be
#' called directly, but to avoid any re-entrant behavior this function is called
#' instead. It is a light wrapper around cogaps_cpp that handles setting
#' the distributed parameters
#' @keywords internal
#'
#' @param index index for which subset to run on
#' @param sets list of all subsets
#' @param data data in a supported format
#' @param allParams list of all parameters
#' @param uncertainty uncertainty of data in the same format
#' @param geneNames names of all genes
#' @param sampleNames names of all samples
#' @param fixedMatrix matrix of matched patterns
#' @return CogapsResult object
callInternalCoGAPS <- function(data, allParams, uncertainty, subsetIndices,
workerID)
{
    # identify which mode of parallelization
    genomeWide <- allParams$gaps@distributed == "genome-wide"
    allParams$gaps@distributed <- NULL

    # subset gene/sample names
    if (genomeWide)
        allParams$geneNames <- allParams$geneNames[subsetIndices]
    else
        allParams$sampleNames <- allParams$sampleNames[subsetIndices]

    allParams$gaps@subsetIndices <- subsetIndices
    allParams$gaps@subsetDim <- ifelse(genomeWide, 1, 2)
    allParams$workerID <- workerID

    # call CoGAPS
    internal <- ifelse(is(data, "character"), cogaps_cpp_from_file, cogaps_cpp)
    raw <- internal(data, allParams, uncertainty)
    return(createCogapsResult(raw, allParams))
}

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
    # randomly sample either rows or columns into subsets to break the data up
    set.seed(allParams$gaps@seed)
    sets <- createSets(data, allParams)
    if (min(sapply(sets, length)) < allParams$gaps@nPatterns)
        stop("data subset dimension less than nPatterns")

    if (is.null(allParams$BPPARAM))
        allParams$BPPARAM <- BiocParallel::MulticoreParam(workers=length(sets))
    
    initialResult <- NULL
    if (is.null(allParams$gaps@fixedPatterns))
    {
        # run Cogaps normally on each subset of the data
        gapsCat(allParams, "Running Across Subsets...\n\n")
        initialResult <- bplapply(1:length(sets), BPPARAM=allParams$BPPARAM,
        FUN=function(i)
        {
            callInternalCoGAPS(data, allParams, uncertainty, sets[[i]], i)
        })

        # get all unmatched patterns
        if (allParams$gaps@distributed == "genome-wide")
            unmatchedPatterns <- lapply(initialResult, function(x) x@sampleFactors)
        else
            unmatchedPatterns <- lapply(initialResult, function(x) x@featureLoadings)

        # match patterns in either A or P matrix
        gapsCat(allParams, "\nMatching Patterns Across Subsets...\n")
        matchedPatterns <- findConsensusMatrix(unmatchedPatterns, allParams$gaps)
    }
    else
    {
        matchedPatterns <- list(consensus=allParams$gaps@fixedPatterns)
    }

    # set fixed matrix
    allParams$gaps@nPatterns <- ncol(matchedPatterns$consensus)
    allParams$gaps@fixedPatterns <- matchedPatterns$consensus
    allParams$gaps@whichMatrixFixed <- ifelse(allParams$gaps@distributed
        == "genome-wide", "P", "A")
        
    # run final phase with fixed matrix
    gapsCat(allParams, "Running Final Stage...\n\n")
    finalResult <- bplapply(1:length(sets), BPPARAM=allParams$BPPARAM,
    FUN=function(i)
    {
        callInternalCoGAPS(data, allParams, uncertainty, sets[[i]], i)
    })

    # concatenate final result
    fullResult <- stitchTogether(finalResult, allParams, sets)

    # add diagnostic information about initial run before returning
    if (!is.null(initialResult)) # check that initial phase was run
    {
        fullResult$diagnostics$firstPass <- initialResult
        fullResult$diagnostics$unmatchedPatterns <- unmatchedPatterns
        fullResult$diagnostics$clusteredPatterns <- matchedPatterns$clusteredPatterns
        fullResult$diagnostics$CorrToMeanPattern <- lapply(matchedPatterns$clusteredPatterns, corrToMeanPattern)
    }

    # include the subsets used
    if (allParams$gaps@distributed == "genome-wide")
        fullResult$diagnostics$subsets <- lapply(sets, function(s) fullResult$geneNames[s])
    else
        fullResult$diagnostics$subsets <- lapply(sets, function(s) fullResult$sampleNames[s])

    # return list, calling function will process this into a CogapsResult object
    return(fullResult)
}


#' find the consensus pattern matrix across all subsets
#' @export
#'
#' @param unmatchedPatterns list of all unmatched pattern matrices from initial
#' run of CoGAPS
#' @param gapsParams list of all CoGAPS parameters
#' @return matrix of consensus patterns
findConsensusMatrix <- function(unmatchedPatterns, gapsParams)
{
    allPatterns <- do.call(cbind, unmatchedPatterns)
    comb <- expand.grid(1:gapsParams@nSets, 1:gapsParams@nPatterns)
    colnames(allPatterns) <- paste(comb[,1], comb[,2], sep=".")
    return(patternMatch(allPatterns, gapsParams))
}

#' Match Patterns Across Multiple Runs
#' @keywords internal
#'
#' @param allPatterns matrix of patterns stored in the columns
#' @param gapsParams CoGAPS parameters object
#' @return a matrix of consensus patterns
#' @importFrom stats weighted.mean
patternMatch <- function(allPatterns, gapsParams)
{
    # cluster patterns
    clusters <- corcut(allPatterns, gapsParams@cut, gapsParams@minNS)

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
    tooLarge <- function(x) ncol(x) > gapsParams@maxNS
    indx <- which(sapply(clusters, tooLarge))
    while (length(indx) > 0)
    {
        clusters <- splitCluster(clusters, indx[1], gapsParams@minNS)
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
#' @return correlation of each pattern
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

    if (any(is.na(corr.dist)))
    {
        stop("NA values in correlation of patterns")
    }

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
#' @param sets indices of sets used to break apart data
#' @return list with all CoGAPS output
stitchTogether <- function(result, allParams, sets)
{
    setIndices <- unlist(sets)
    if (allParams$gaps@distributed == "genome-wide")
    {
        # combine A matrices, re-order so it matches original data
        Amean <- do.call(rbind, lapply(result, function(x) x@featureLoadings))
        Asd <- do.call(rbind, lapply(result, function(x) x@featureStdDev))

        # copy P matrix - same for all sets
        Pmean <- result[[1]]@sampleFactors
        Psd <- matrix(0, nrow=nrow(Pmean), ncol=ncol(Pmean))

        # if each feature was used once, re-order to match data
        if (nrow(Amean) == length(setIndices))
        {
            indices <- 1:nrow(Amean)
            if (identical(sort(indices), sort(setIndices)))
            {
                reorder <- match(indices, setIndices)
                Amean <- Amean[reorder,]
                Asd <- Asd[reorder,]
            }
        }
    }
    else
    {
        # combine P matrices, re-order so it matches original data
        Pmean <- do.call(rbind, lapply(result, function(x) x@sampleFactors))
        Psd <- do.call(rbind, lapply(result, function(x) x@sampleStdDev))

        # copy A matrix - same for all sets
        Amean <- result[[1]]@featureLoadings
        Asd <- matrix(0, nrow=nrow(Amean), ncol=ncol(Amean))

        # if each sample was used once, re-order to match data
        if (nrow(Pmean) == length(setIndices))
        {
            indices <- 1:nrow(Pmean)
            if (identical(sort(indices), sort(setIndices)))
            {
                reorder <- match(indices, setIndices)
                Pmean <- Pmean[reorder,]
                Psd <- Psd[reorder,]
            }
        }
    }

    return(list("Amean"=Amean, "Asd"=Asd, "Pmean"=Pmean, "Psd"=Psd,
        "seed"=allParams$gaps@seed, "geneNames"=rownames(Amean),
        "sampleNames"=rownames(Pmean),
        "meanChiSq"=sum(sapply(result, function(r) r@metadata$meanChiSq))))
}

