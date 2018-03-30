
#' patternMatch4Parallel
#'
#' @param Atot a matrix containing the total by set estimates of Pmean output from \code{reOrderBySet}
#' @param nSets number of parallel sets used to generate \code{Atot}
#' @param cnt  number of branches at which to cut dendrogram
#' @param minNS minimum of individual set contributions a cluster must contain
#' @param ignore.NA logical indicating whether or not to ignore NAs from potential over dimensionalization. Default is FALSE.
#' @param bySet logical indicating whether to return list of matched set solutions from \code{Atot}
#' @param R2mean 
#' @param
#' @param ... additional parameters for \code{agnes}
#' @return a matrix of concensus patterns by samples. If \code{bySet=TRUE} then a list of the set contributions to each
#' concensus pattern is also returned.
#' @export
#' @seealso \code{\link{fastcluster}}
#'
#'

cellMatchR <- function(Atot,nSets, cnt, minNS=NULL, maxNS=NULL, ignore.NA=FALSE, bySet=FALSE, plotDen=FALSE,...){

  if(is.null(minNS)){minNS=nSets/2}
  if(is.null(maxNS)){maxNS=nSets+minNS}

  if(ignore.NA==FALSE){if(anyNA(Atot)){
    warning('Non-sparse matrixes produced. Reducing the number of patterns asked for and rerun.')
  }}
  if(ignore.NA==TRUE){Atot<-Atot[complete.cases(Atot),]}


corcut<-function(Atot,minNS,cnt,cluster.method){
  corr.dist=cor(Atot)
  corr.dist=1-corr.dist

  clust=agnes(x=corr.dist,diss=TRUE,cluster.method)
  #clust=fastcluster::hclust(dist(corr.dist))
  cut=cutree(as.hclust(clust),k=cnt)

  cls=sort(unique(cut))
  cMNs=matrix(ncol=cnt,nrow=dim(Atot)[1])
  colnames(cMNs)=cls
  rownames(cMNs)=rownames(Atot)

  RtoMeanPattern <- list()
  AByClust <- list()
    for(i in cls){
        if (is.null(dim(Atot[,cut == i]))==TRUE){
          next
        } else if(dim(Atot[,cut == i])[2] < minNS){
            next
        } else{
            cMNs[,i]=rowMeans(Atot[,cut==i])
            AByClust[[i]] <- Atot[,cut==i]
            nIN=sum(cut==i)
            RtoMeanPattern[[i]] <- sapply(1:nIN,function(j) {round(cor(x=Atot[,cut==i][,j],y=cMNs[,i]),3)})
        }
      }
    PByClust[sapply(PByClust,is.null)]<-NULL
    RtoMeanPattern[sapply(RtoMeanPattern,is.null)]<-NULL
    return(list("RtoMeanPattern"=RtoMeanPattern,"AByClust"=AByClust))
  }   

  cc<-corcut(Atot,minNS,cnt,cluster.method)

    ### split by maxNS
    indx<-which(unlist(lapply(cc$PByClust,function(x) dim(x)[1]>maxNS)))
    while(length(indx)>0){
          icc<-corcut(cc$AByClust[[indx[1]]],minNS,2,cluster.method)
          cc$AByClust[[indx[1]]]<-icc[[2]][[2]]
          cc$RtoMeanPattern[[indx[1]]]<-icc[[1]][[2]]
          if(length(icc[[2]])>1){
                cc$AByClust<-append(cc$AByClust,icc[[2]][1])
                cc$RtoMeanPattern<-append(cc$RtoMeanPattern,icc[[1]][1])
          }
      indx<-which(unlist(lapply(cc$PByClust,function(x) dim(x)[1]>maxNS)))
    }

#weighted.mean(AByClustDrop[[1]],RtoMPDrop[[1]])
AByCDSWavg<- t(sapply(1:length(cc$AByClust),function(z) apply(cc$AByClust[[z]],1,function(x) weighted.mean(x,(cc$RtoMeanPattern[[z]])^3))))
rownames(AByCDSWavg) <- lapply(1:length(cc$AByClust),function(x) paste("Pattern",x))

#scale As
Amax <- apply(AByCDSWavg,1,max)
AByCDSWavgScaled <- t(sapply(1:dim(AByCDSWavg)[1],function(x) AByCDSWavg[x,]/Amax[x]))
rownames(AByCDSWavgScaled) <- rownames(AByCDSWavg)

  if(bySet){  
    return(list("consenusAs"=t(AByCDSWavgScaled),"ABySet"=cc))
  } else {return(AByCDSWavgScaled)}

}
