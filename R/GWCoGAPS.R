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
#'@param fname character string used to label file output. Default is "GWCoGAPS.out"
#'@paramPatternsMatchFN function to use for pattern matching across sets
#'@paramCut number of branches at which to cut dendrogram used in patternMatch4Parallel
#'@paramminNS minimum of individual set contributions a cluster must contain
#'@param ABins a matrix of same size as A which gives relative
#' probability of that element being non-zero
#'@param PBins a matrix of same size as P which gives relative
#' probability of that element being non-zero
#'@param  simulation_id name to attach to atoms files if created
#'@param  nEquil number of iterations for burn-in
#'@param  nSample number of iterations for sampling
#'@param  nOutR how often to print status into R by iterations
#'@param  output_atomic whether to write atom files (large)
#'@param sampleSnapshots Boolean to indicate whether to capture
#' individual samples from Markov chain during sampling
#'@param numSnapshots the number of individual samples to capture
#'@param alphaA sparsity parameter for A domain
#'@param max_gibbmass_paraA limit truncated normal to max size
#'@param alphaP sparsity parameter for P domain
#'@param max_gibbmass_paraP limit truncated normal to max size
#'@return list with A and P matrix estimates, chi-squared and atom
#'         numbers of sample by iteration, and chi-squared of mean
#'@export
#'@seealso \code{\link{CoGAPS}}
#' @examples \dontrun{
#' GWCoGAPS(nCores=NA, D, S, nFactor, nSets,saveBySetResults=TRUE, fname=fname,
#' PatternsMatchFN = patternMatch4Parallel,numSnapshots=numSnapshots,minNS=minNS)
#' }
#'

GWCoGAPS <- function(D, S, nFactor, nSets, nCores=NA,
                                 saveBySetResults=FALSE, fname="GWCoGAPS.out",
                                 PatternsMatchFN = postCoGAPSPatternMatch,
                                 manualMatch=FALSE, Cut=NA, minNS=NA,
                                 numSnapshots=numSnapshots, ...) {

# establish the number of cores that you are able to use
if(is.na(nCores)){nCores<-nSets}
registerDoParallel(cores=nCores)

# break the data into sets
genesInSets<-createGWCoGAPSSets(data=D, nSet=nSets, keep=FALSE)

# set gene min to 0
D <- sweep(D, 1, apply(D,1,function(x){pmin(0,min(x))}))

#generate seeds for parallelization
nut<-generateSeeds(chains=nSets, seed=-1)

# run CoGAPS for each set
AP <- foreach(i=1:nSets) %dopar% {
     D <- as.matrix(D[genesInSets[[i]],])
     S <- as.matrix(S[genesInSets[[i]],])
     gapsRun(D=D, S=S, nFactor=nFactor,seed=nut[i],numSnapshots=numSnapshots)
}

BySet<-reOrderBySet(AP=AP,nFactor=nFactor,nSets=nSets)

#run postpattern match function
if(is.na(Cut)){Cut<-nFactor}
matchedPs<-patternMatch4Parallel(Ptot=BySet$P,nP=nFactor,nSets=nSets,cnt=Cut,minNS=minNS,bySet=TRUE)

#save BySet outputs
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
	AP <- gapsMapRun(D, S, FP=matchedPs, nFactor=nFactorFinal, fixedMatrix = "P",seed=nut[i],numSnapshots=numSnapshots)
    }

#extract A and Asds
As4fixPs <- postFixed4Parallel(AP.fixed=Fixed,setPs=matchedPs)


#save final
AP.fixed<-list("A"=As4fixPs$A, "Asd"=As4fixPs$Asd, "P"=matchedPs,"PbySet"=PbySet)
save(AP.fixed, file=paste("WGCoGAPS.AP.fixed.",fname,".Rda",sep=""))
print(paste("WGCoGAPS.AP.fixed.",fname,".Rda",sep=""))
}
