#' patternMatch4Parallel
#'
#' @param Ptot a matrix containing the total by set estimates of Pmean output
#' from \code{reOrderBySet}
#' @param nSets number of parallel sets used to generate \code{Ptot}
#' @param cnt  number of branches at which to cut dendrogram
#' @param minNS minimum of individual set contributions a cluster must contain
#' @param maxNS max of individual set contributions a cluster must contain.
#' default is nSets+minNS
#' @param cluster.method the agglomeration method to be used for clustering
#' @param ignore.NA logical indicating whether or not to ignore NAs from
#' potential over dimensionalization. Default is FALSE.
#' @param bySet logical indicating whether to return list of matched set
#' solutions from \code{Ptot}
#' @param ... additional parameters for \code{agnes}
#' @return a matrix of concensus patterns by samples. If \code{bySet=TRUE} then
#' a list of the set contributions to each
#' concensus pattern is also returned.
#' @seealso \code{\link{agnes}}
#' @export
patternMatch4Parallel <- function(Ptot, nSets, cnt, minNS=NA, maxNS=NULL, 
cluster.method="complete", ignore.NA=FALSE, bySet=FALSE, ...)
{
    if (is.na(minNS))
    {
        minNS <- ceiling(nSets / 2)
    }
    if (is.null(maxNS))
    {
        maxNS=nSets+minNS
    }
    if (ignore.NA==FALSE & anyNA(Ptot))
    {
        warning(paste("Non-sparse matrixes produced. Reducing the number of",
            "patterns asked for and rerun."))
    }
    if (ignore.NA==TRUE)
    {
        Ptot <- Ptot[complete.cases(Ptot),]
    }

    corcut <-function(Ptot,minNS,cnt,cluster.method)
    {
        # corr dist
        corr.dist <- cor(t(Ptot))
        corr.dist <- 1-corr.dist
        # cluster
        #library(cluster)
        clust <- agnes(x=corr.dist, diss=TRUE, method=cluster.method)
        cut <- cutree(as.hclust(clust), k=cnt)
        #save.image(file=paste("CoGAPS.", nP, "P.", nS, "Set.CorrClustCut",
            #cnt,".RData"))

        cls <- sort(unique(cut))
        cMNs <- matrix(nrow=cnt, ncol=dim(Ptot)[2])
        rownames(cMNs) <- cls
        colnames(cMNs) <- colnames(Ptot)

        RtoMeanPattern <- list()
        PByClust <- list()
        for (i in cls)
        {
            if (!is.null(dim(Ptot[cut == i,])))
            {
                if (dim(Ptot[cut == i,])[1] >= minNS)
                {
                    cMNs[i,] <- colMeans(Ptot[cut==i,])
                    PByClust[[i]] <- Ptot[cut==i,]
                    nIN <- sum(cut==i)
                    RtoMeanPattern[[i]] <- sapply(1:nIN, function(j)
                        round(cor(x=Ptot[cut==i,][j,], y=cMNs[i,]),3)
                    )
                }
            }
        }

        PByClust[sapply(PByClust,is.null)] <- NULL
        RtoMeanPattern[sapply(RtoMeanPattern,is.null)] <- NULL
        return(list("RtoMeanPattern"=RtoMeanPattern, "PByClust"=PByClust))
    }    

    cc <- corcut(Ptot,minNS,cnt,cluster.method)

    ### split by maxNS
    indx <- which(unlist(lapply(cc$PByClust,function(x) dim(x)[1]>maxNS)))
    while(length(indx) > 0)
    {
        icc <- corcut(cc$PByClust[[indx[1]]],minNS,2,cluster.method)
        if (length(icc[[2]])==0)
        {
          indx <- indx[-1]
          next
        }
        else
        {
            cc$PByClust[[indx[1]]] <- icc[[2]][[1]]
            cc$RtoMeanPattern[[indx[1]]] <- icc[[1]][[1]]
            if (length(icc[[2]]) > 1)
            {
                cc$PByClust <- append(cc$PByClust,icc[[2]][2])
                cc$RtoMeanPattern <- append(cc$RtoMeanPattern,icc[[1]][2])
            } 
            indx <- which(unlist(lapply(cc$PByClust, function(x)
                dim(x)[1] > maxNS)))
        }
    }

    #weighted.mean
    PByCDSWavg <- t(sapply(1:length(cc$PByClust), function(z)
        apply(cc$PByClust[[z]],2, function(x)
            weighted.mean(x, (cc$RtoMeanPattern[[z]])^3))))
    rownames(PByCDSWavg) <- lapply(1:length(cc$PByClust), function(x)
        paste("Pattern", x))

    #scale ps
    Pmax <- apply(PByCDSWavg,1,max)
    PByCDSWavgScaled <- t(sapply(1:dim(PByCDSWavg)[1], function(x)
        PByCDSWavg[x,] / Pmax[x]))
    rownames(PByCDSWavgScaled) <- rownames(PByCDSWavg)

    if (bySet)
    {
        return(list("consenusPatterns"=PByCDSWavgScaled, "PBySet"=cc))
    }
    else
    {
        return("consenusPatterns"=PByCDSWavgScaled)
    }
}
