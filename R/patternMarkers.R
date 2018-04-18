#' patternMarkers
#'
#' @param Amatrix A matrix of genes by weights resulting from CoGAPS or other NMF decomposition
#' @param scaledPmatrix logical indicating whether the corresponding pattern matrix was fixed to have max 1 during decomposition
#' @param Pmatrix the corresponding Pmatrix (patterns X samples) for the provided Amatrix (genes x patterns). This must be supplied if scaledPmatrix is FALSE.
#' @param threshold # the type of threshold to be used. The default "all" will distribute genes into pattern with the lowest ranking. The "cut" thresholding by the first gene to have a lower ranking, i.e. better fit to, a pattern.
#' @param lp a vector of weights for each pattern to be used for finding markers. If NA markers for each pattern of the A matrix will be used.
#' @param full logical indicating whether to return the ranks of each gene for each pattern
#' @return By default a non-overlapping list of genes associated with each \code{lp}. If \code{full=TRUE} a data.frame of
#' genes rankings with a column for each \code{lp} will also be returned.
#' @export
patternMarkers <- function(Amatrix=NA, scaledPmatrix=FALSE, Pmatrix=NA,
threshold="all", lp=NA, full=FALSE)
{
    if(scaledPmatrix==FALSE)
    {
        if(!is.na(Pmatrix)){
        pscale <- apply(Pmatrix,1,max)   # rescale p's to have max 1
        Amatrix <- sweep(Amatrix, 2, pscale, FUN="*")   # rescale A in accordance with p's having max 1
        }
        else(warning("P values must be provided if not already scaled"))
    }

    # find the A with the highest magnitude
    Arowmax <- t(apply(Amatrix, 1, function(x) x/max(x)))
        
    # determine which genes are most associated with each pattern
    sstat<-matrix(NA, nrow=nrow(Amatrix), ncol=ncol(Amatrix),dimnames=dimnames(Amatrix))
    ssranks<-matrix(NA, nrow=nrow(Amatrix), ncol=ncol(Amatrix),dimnames=dimnames(Amatrix))#list()
    ssgenes<-matrix(NA, nrow=nrow(Amatrix), ncol=ncol(Amatrix),dimnames=NULL)
    nP=dim(Amatrix)[2]
    if(!is.na(lp))
    {
        if(length(lp)!=dim(Amatrix)[2]){
            warning("lp length must equal the number of columns of the Amatrix")
        }
        for (i in 1:nP){
            sstat[,i] <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
            ssranks[,i]<-rank(sstat[,i])
            ssgenes[,i]<-names(sort(sstat[,i],decreasing=FALSE,na.last=TRUE))
        }
    }
    else
    {
        for(i in 1:nP){
            lp <- rep(0,dim(Amatrix)[2])
            lp[i] <- 1
            sstat[,i] <- unlist(apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp))))
            #ssranks[order(sstat[,i]),i] <- 1:dim(sstat)[1]
            ssranks[,i]<-rank(sstat[,i])
            ssgenes[,i]<-names(sort(sstat[,i],decreasing=FALSE,na.last=TRUE))
        }
    }

    if(threshold=="cut"){
        geneThresh <- sapply(1:nP,function(x) min(which(ssranks[ssgenes[,x],x] > apply(ssranks[ssgenes[,x],],1,min))))
        ssgenes.th <- sapply(1:nP,function(x) ssgenes[1:geneThresh[x],x])
    }
    else if (threshold=="all"){
        pIndx<-unlist(apply(sstat,1,which.min))
        gBYp <- list()
        for(i in sort(unique(pIndx))){
            gBYp[[i]]<-sapply(strsplit(names(pIndx[pIndx==i]),"[.]"),function(x) x[[1]][1])

        }
        ssgenes.th <- lapply(1:max(sort(unique(pIndx))), function(x) {
            ssgenes[which(ssgenes[,x] %in% gBYp[[x]]),x]
        })
    }
    else{
        stop("Threshold arguement not viable option")
    }

    if (full)
        return(list("PatternMarkers"=ssgenes.th,"PatternRanks"=ssranks,"PatternMarkerScores"=sstat))
    else
        return("PatternMarkers"=ssgenes.th)
}
