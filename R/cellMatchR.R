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
#' @param plotDen plot
#' @param ... additional parameters for \code{agnes}
#' @return a matrix of concensus patterns by samples. If \code{bySet=TRUE} then a list of the set contributions to each
#' concensus pattern is also returned.
cellMatchR <- function(Atot, nSets, cnt, minNS=NA, maxNS=NA,
ignore.NA=FALSE, bySet=FALSE, plotDen=FALSE, ...)
{
    if (is.na(minNS))
    {
        minNS <- nSets / 2
    }
    if (is.na(maxNS))
    {
        maxNS <- nSets + minNS
    }
    if (!ignore.NA)
    {
        if (anyNA(Atot))
        {
            warning(paste("Non-sparse matrixes produced. Reducing the number",
                "of patterns asked for and rerun."))
        }
    }
    else
    {
        Atot <- Atot[complete.cases(Atot),]
    }

    corcut <- function(Atot, minNS, cnt, cluster.method)
    {
        corr.dist <- cor(Atot)
        corr.dist <- 1-corr.dist

        clust <- agnes(x=corr.dist, diss=TRUE, cluster.method)
        #clust=fastcluster::hclust(dist(corr.dist))
        cut <- cutree(as.hclust(clust),k=cnt)

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

        AByClust[sapply(AByClust,is.null)] <- NULL
        RtoMeanPattern[sapply(RtoMeanPattern,is.null)] <- NULL
        return(list("RtoMeanPattern"=RtoMeanPattern, "AByClust"=AByClust))
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

    if (bySet)
    {  
        return(list("consenusAs"=t(AByCDSWavgScaled),"ABySet"=cc))
    }
    else
    {
        return(AByCDSWavgScaled)
    }
}

# pattern first then file name
# dont overwrite
# some elements are negative