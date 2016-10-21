
#' patternMatch4Parallel
#'
#' @param Ptot a matrix containing the total by set estimates of Pmean output from \code{reOrderBySet}
#' @param nSets number of parallel sets used to generate \code{Ptot}
#' @param cnt  number of branches at which to cut dendrogram
#' @param minNS minimum of individual set contributions a cluster must contain
#' @param cluster.method the agglomeration method to be used for clustering
#' @param ignore.NA logical indicating whether or not to ignore NAs from potential over dimensionalization. Default is FALSE.
#' @param bySet logical indicating whether to return list of matched set solutions from \code{Ptot}
#' @param ... additional parameters for \code{agnes}
#'
#' @return a matrix of concensus patterns by samples. If \code{bySet=TRUE} then a list of the set contributions to each
#' concensus pattern is also returned.
#' @export
#' @seealso \code{\link{agnes}}
#'
#'

patternMatch4Parallel <- function(Ptot,
  nSets, #number of sets
  cnt, # number of branches at which to cut dendrogram
  minNS, # minimum of sets a cluster must contain
  cluster.method="complete",
  ignore.NA=FALSE, # ignore NAs from potential over dimensionalization
  bySet=FALSE,
  ...){

#### read in CoGAPS results
cdir <- getwd()
#if(!is.null(path)){setwd(path)}
if(!is.null(minNS)){minNS=nSets/2}

if(ignore.NA==FALSE){if(anyNA(Ptot)){warning('Non-sparse matrixes produced. Reducing the number of patterns asked for and rerun.')}}
if(ignore.NA==TRUE){Ptot<-Ptot[complete.cases(Ptot),]}

####################################################################
# corr dist
corr.dist=cor(t(Ptot))
corr.dist=1-corr.dist
# cluster
clust=agnes(x=corr.dist,diss=T,method=cluster.method)
cut=cutree(as.hclust(clust),k=cnt)
#save.image(file=paste("CoGAPS.",nP,"P.",nS,"Set.CorrClustCut",cnt,".RData"))
####################################################################
#drop n<4 and get weighted Avg
cls=sort(unique(cut))
cMNs=matrix(nrow=cnt,ncol=dim(Ptot)[2])
rownames(cMNs)=cls
colnames(cMNs)=colnames(Ptot)

RtoMeanPattern <- list()
PByClust <- list()
for(i in cls){
   if(is.null(dim(Ptot[cut == i, ]))==TRUE){
       cMNs[i,] <- Ptot[cut == i, ]
       RtoMeanPattern[[i]] <- rep(1,length(Ptot[cut == i, ]))
       PByClust[[i]] <- t(as.matrix(Ptot[cut == i, ]))
   }
  else{
  cMNs[i,]=colMeans(Ptot[cut==i,])
  PByClust[[i]] <- Ptot[cut==i,]
  nIN=sum(cut==i)
  RtoMeanPattern[[i]] <- sapply(1:nIN,function(j) {round(cor(x=Ptot[cut==i,][j,],y=cMNs[i,]),3)})
  }
}

#drop n<4 and drop less than .7
PByClustDrop <- list()
RtoMPDrop <- list()
for(i in cls){
  if(is.null(dim(PByClust[[i]]))==TRUE){next}
  if(dim(PByClust[[i]])[1]<minNS){next}
  else{
    indx <- which(RtoMeanPattern[[i]]>.7,arr.ind = TRUE)
    PByClustDrop <- append(PByClustDrop,list(PByClust[[i]][indx,]))
    RtoMPDrop <- append(RtoMPDrop,list(RtoMeanPattern[[i]][indx]))
  }
}


### split by corr  (build in drop if below 4)
PByCDS <- list()
RtoMPDS <- list()
for(j in 1:length(PByClustDrop)){
  if(is.null(dim(PByClustDrop[[j]]))==TRUE){
      next
      }
  if(dim(PByClustDrop[[j]])[1]<minNS+nSets){
    PByCDS <- append(PByCDS,PByClustDrop[j])
    RtoMPDS <- append(RtoMPDS,RtoMPDrop[j])
  }
  if(dim(PByClustDrop[[j]])[1]>=minNS+nSets){
    corr.distPBCD=cor(t(PByClustDrop[[j]]))
    corr.distPBCD=1-corr.distPBCD
    library(cluster)
    clustPBCD=agnes(x=corr.distPBCD,diss=T,method="complete")
    cutPBCD=cutree(as.hclust(clustPBCD),k=2)
    g1 <- PByClustDrop[[j]][cutPBCD==1,]
    PByCDS <- append(PByCDS,list(g1))
    RtoMPDS <- append(RtoMPDS,list(sapply(1:dim(g1)[1],function(z) round(cor(x=g1[z,],y=colMeans(PByClustDrop[[j]][cutPBCD==1,])),3))))
    g2 <- PByClustDrop[[j]][cutPBCD==2,]
    if (is.null(dim(g2)[1])==FALSE){
        PByCDS <- append(PByCDS,list(g2))
        RtoMPDS <- append(RtoMPDS,list(sapply(1:dim(g2)[1],function(z) round(cor(x=g2[z,],y=colMeans(PByClustDrop[[j]][cutPBCD==2,])),3))))
    }
  }
#print(j)
#print(str(PByCDS))
}

#weighted.mean(PByClustDrop[[1]],RtoMPDrop[[1]])
PByCDSWavg<- t(sapply(1:length(PByCDS),function(z) apply(PByCDS[[z]],2,function(x) weighted.mean(x,(RtoMPDS[[z]])^3))))
rownames(PByCDSWavg) <- lapply(1:length(PByCDS),function(x) paste("Pattern",x))

#save
#save(PByCDSWavg,file=paste("PAatternSummary.UnScaled.",fname,".CoGAPS.",nP,"P.",nS,"Set.CorrClustCut",cnt,".RData",sep=""))

#scale ps
Pmax <- apply(PByCDSWavg,1,max)
PByCDSWavgScaled <- t(sapply(1:dim(PByCDSWavg)[1],function(x) PByCDSWavg[x,]/Pmax[x]))
rownames(PByCDSWavgScaled) <- rownames(PByCDSWavg)

if(bySet){
# return by set and final
PBySet<-PByCDS
return(list("consenusPatterns"=PByCDSWavgScaled,"PBySet"=PBySet))
} else {return(PByCDSWavgScaled)}

}


####################################################################
####################################################################

