#' GWCoGAPS
#'
#'\code{GWCoGAPS} calls the C++ MCMC code and performs Bayesian
#'matrix factorization returning the two matrices that reconstruct
#'the data matrix for whole genome data;
#'
#'@param D data matrix
#'@param S uncertainty matrix (std devs for chi-squared of Log Likelihood)
#'@param nFactor number of patterns (basis vectors, metagenes), which must be
#'greater than or equal to the number of rows of FP
#'@param nSets number of sets for parallelization
#'@param nCores number of cores for parallelization. If left to the default NA, nCores = nSets.
#'@param saveBySetResults logical indicating whether to save by intermediary by set results. Default is FALSE.
#'@param fname character string used to label file output. Default is "GWCoGAPS.AP.fixed"
#'@param PatternsMatchFN function to use for pattern matching across sets
#'@param Cut number of branches at which to cut dendrogram used in patternMatch4Parallel
#'@param minNS minimum of individual set contributions a cluster must contain
#'@param ... additional parameters to be fed into \code{gapsRun} and \code{gapsMapRun}
#'@export
#'@seealso \code{\link{gapsRun}}, \code{\link{patternMatch4Parallel}}, and \code{\link{gapsMapRun}}
#' @examples \dontrun{
#' GWCoGAPS(nCores=NA, D, S, nFactor, nSets,saveBySetResults=TRUE, fname=fname,
#' PatternsMatchFN = patternMatch4Parallel,numSnapshots=numSnapshots,minNS=minNS)
#' }
#'

GWCoGAPS <- function(D, S, nFactor, nSets, nCores=NA,
                                 saveBySetResults=FALSE, fname="GWCoGAPS.AP.fixed",
                                 PatternsMatchFN = patternMatch4Parallel,
                                 Cut=NA, minNS=NA, ...) {

# establish the number of cores that you are able to use
if(is.na(nCores)){nCores<-nSets}
registerDoParallel(cores=nCores)

# break the data into sets
genesInSets<-createGWCoGAPSSets(data=D, nSets=nSets, keep=FALSE)

# set gene min to 0
D <- sweep(D, 1, apply(D,1,function(x){pmin(0,min(x))}))

#generate seeds for parallelization
nut<-generateSeeds(chains=nSets, seed=-1)

# run CoGAPS for each set
AP <- foreach(i=1:nSets) %dopar% {
     D <- as.matrix(D[genesInSets[[i]],])
     S <- as.matrix(S[genesInSets[[i]],])
     gapsRun(D=D, S=S, nFactor=nFactor,seed=nut[i],...)
}

BySet<-reOrderBySet(AP=AP,nFactor=nFactor,nSets=nSets)

#run postpattern match function
if(is.na(Cut)){Cut<-nFactor}
matchedPs<-patternMatch4Parallel(Ptot=BySet$P,nP=nFactor,nSets=nSets,cnt=Cut,minNS=minNS,bySet=TRUE)

#save BySet outputs
class(AP)<-append(class(AP),"CoGAPS")
if(saveBySetResults==TRUE){
    save(AP,BySet,matchedPs,file=sprintf('APbySet.%s.nP%d.set%d.Rda',fname, nFactor,nSets))
    print(sprintf('APbySet.%s.nP%d.set%d.Rda',fname, nFactor,nSets))
}

PbySet<-matchedPs[["PBySet"]]
matchedPs<-matchedPs[[1]]

#generate seeds for parallelization
nut<-generateSeeds(chains=nSets, seed=-1)

#final number of factors
nFactorFinal<-dim(matchedPs)[1]

# run fixed CoGAPS
Fixed <- foreach(i=1:nSets) %dopar% {
	D <- as.matrix(D[genesInSets[[i]],])
	S <- as.matrix(S[genesInSets[[i]],])
	AP <- gapsMapRun(D, S, FP=matchedPs, nFactor=nFactorFinal, fixedMatrix = "P",seed=nut[i],...)
    }

#extract A and Asds
As4fixPs <- postFixed4Parallel(AP.fixed=Fixed,setPs=matchedPs)


#save final
AP.fixed<-list("Amean"=As4fixPs$A, "Asd"=As4fixPs$Asd, "Pmean"=matchedPs,"PbySet"=PbySet)
class(AP.fixed)<-append(class(AP.fixed),"CoGAPS")
save(AP.fixed, file=paste(fname,".Rda",sep=""))
print(paste(fname,".Rda",sep=""))
}
