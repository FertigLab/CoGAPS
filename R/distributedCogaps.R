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
    initialA <- initialP <- matrix(0)
    if (allParams$modelParams@distributed == "genome-wide")
    {
        allParams$whichMatrixFixed <- "P"
        initialP <- consensusMatrix
    }
    else
    {
        allParams$whichMatrixFixed <- "A"
        initialA <- consensusMatrix
    }

    # run with a fixed matrix
    finalResult <- foreach(i=1:nSets) %dopar%
    {
        cogaps_cpp(data, allParams, uncertainty, sets[[i]], initialA, initialP)
    }

    # get result 
    resultList <- stitchTogether(finalResult, allParams)
    resultList$seed <- allParams$modelParams@seed
    resultList$meanChiSq <- sum(sapply(finalResult, function(r) r$meanChiSq))
    resultList$diagnostics <- list()
    return(resultList)
}

nrow_helper <- function(matrix)
{
    if (class(matrix) == "character")
    {
        switch(file_ext(matrix),
            "csv" = nrow(data.table::fread(matrix, select=1)),
            "tsv" = nrow(data.table::fread(matrix, select=1)),
            "mtx" = as.numeric(data.table::fread(matrix, nrows=1, fill=TRUE)[1,1])
        )
    }
    return(nrow(matrix))
}

ncol_helper <- function(matrix)
{
    if (class(matrix) == "character")
    {
        switch(file_ext(matrix),
            "csv" = ncol(data.table::fread(matrix, nrows=1)) - 1,
            "tsv" = ncol(data.table::fread(matrix, nrows=1)) - 1,
            "mtx" = as.numeric(data.table::fread(matrix, nrows=1, fill=TRUE)[1,2])
        )
    }
    return(ncol(matrix))
}        

createSets <- function(data, allParams)
{
    total <- ifelse(allParams$modelParams@distributed == "genome-wide",
        nrow_helper(data), ncol_helper(data))
    all_indices <- 1:total
    setSize <- floor(total / nSets)
    sets <- list()
    for (n in 1:(allParams$gaps@nSets - 1))
    {
        selected <- sample(all_indices, setSize, replace=FALSE)
        sets[[n]] <- selected
        all_indices <- setdiff(all_indices, selected)
    }
    sets[[allParams$gaps@nSets]] <- all_indices
    return(sets)
}

findConsensusMatrix <- function(result, allParams)
{
    if (allParams$modelParams@distributed == "genome-wide")
        allPatterns <- do.call(cbind, lapply(initialResult, function(x) x@sampleFactors))
    else
        allPatterns <- do.call(cbind, lapply(initialResult, function(x) x@featureLoadings))

    comb <- expand.grid(1:allParams$gaps@nSets, 1:allParams$gaps@nPatterns)
    rownames(allPatterns) <- paste(comb[,1], comb[,2], sep=".")
    matched <- patternMatch(allPatterns, allParams)
}

#' cellMatchR
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
#' @return a matrix of concensus patterns by samples. If \code{bySet=TRUE} then a list of the set contributions to each
#' concensus pattern is also returned.
patternMatch <- function(matrix, allParams)
{
    if (is.na(minNS))
    {
        minNS <- nSets / 2
    }
    if (is.na(maxNS))
    {
        maxNS <- nSets + minNS
    }

    cc <- corcut(Atot, minNS, cnt, cluster.method)

    ### split by maxNS
    indx <- which(unlist(lapply(cc$AByClust, function(x) dim(x)[1] > maxNS)))
    i <- 1
    while (length(indx) > 0)
    { 
        icc <- corcut(cc$AByClust[[indx[1]]], minNS, 2, cluster.method)
        if (length(icc[[2]]) == 0)
        {
            indx <- indx[-1]
            next
        }
        else
        {
            cc$AByClust[[indx[1]]] <- icc[[2]][[1]]
            cc$RtoMeanPattern[[indx[1]]] <- icc[[1]][[1]]
            if (length(icc[[2]]) > 1)
            {
                cc$AByClust<-append(cc$AByClust,icc[[2]][2])
                cc$RtoMeanPattern<-append(cc$RtoMeanPattern,icc[[1]][2])
            } 
            indx <- which(unlist(lapply(cc$AByClust,function(x) dim(x)[1]>maxNS)))
        }
    }

    #weighted.mean(AByClustDrop[[1]],RtoMPDrop[[1]])
    AByCDSWavg <- t(sapply(1:length(cc$AByClust), function(z)
        apply(cc$AByClust[[z]], 1, function(x)
            weighted.mean(x, (cc$RtoMeanPattern[[z]])^3))))
    rownames(AByCDSWavg) <- lapply(1:length(cc$AByClust), function(x)
        paste("Pattern",x))
    #scale As
    Amax <- apply(AByCDSWavg, 1, max)
    AByCDSWavgScaled <- t(sapply(1:dim(AByCDSWavg)[1], function(x)
        AByCDSWavg[x,] / Amax[x]))
    rownames(AByCDSWavgScaled) <- rownames(AByCDSWavg)

    return(AByCDSWavgScaled)
}

corcut <- function(matrix, allParams)
{
    corr.dist <- cor(Atot)
    corr.dist <- 1-corr.dist

    clust <- agnes(x=corr.dist, diss=TRUE, "complete")
    cut <- cutree(as.hclust(clust), k=cnt)

    cls <- sort(unique(cut))
    cMNs <- matrix(ncol=cnt,nrow=dim(Atot)[1])
    colnames(cMNs) <- cls
    rownames(cMNs) <- rownames(Atot)

    RtoMeanPattern <- list()
    AByClust <- list()
    for (i in cls)
    {
        if (!is.null(dim(Atot[,cut == i])))
        {
            if (dim(Atot[,cut == i])[2] >= minNS)
            {
                cMNs[,i] <- rowMeans(Atot[,cut==i])
                AByClust[[i]] <- Atot[,cut==i]
                nIN <- sum(cut==i)
                RtoMeanPattern[[i]] <- sapply(1:nIN, function(j) 
                    round(cor(x=Atot[,cut==i][,j], y=cMNs[,i]),3))
            }
        }
    }
    return(list("RtoMeanPattern"=RtoMeanPattern, "AByClust"=AByClust))
}   

stitchTogether <- function(result, allParams)
{
    if (allParams$modelParams@distributed == "genome-wide")
    {
        return(list(
            "Amean" = do.call(rbind, lapply(AP.fixed, function(x) x@featureLoadings)),
            "Asd"   = do.call(rbind, lapply(AP.fixed, function(x) x@featureStdDev)),
            "Pmean" = consensusPatterns,
            "Psd"   = matrix(0, nrow=nrow(P), ncol=ncol(P))
        ))
    }
    else
    {
        return(list(
            "Amean" = consensusPatterns,
            "Asd"   = matrix(0, nrow=nrow(A), ncol=ncol(A)),
            "Pmean" = do.call(rbind, lapply(AP.fixed, function(x) x@sampleFactors)),
            "Psd"   = do.call(rbind, lapply(AP.fixed, function(x) x@sampleStdDev))
        ))
    }
}

