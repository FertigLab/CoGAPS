#' Genome Wide CoGAPS
#'
#' @details calls the C++ MCMC code and performs Bayesian
#'  matrix factorization returning the two matrices that reconstruct
#'  the data matrix for whole genome data;
#' @param D data matrix
#' @param S uncertainty matrix (std devs for chi-squared of Log Likelihood)
#' @param params GapsParams object 
#' @param nSets number of sets for parallelization
#' @param nCores number of cores for parallelization. If left to the default NA, nCores = nSets.
#' @param saveBySetResults logical indicating whether to save by intermediary by set results. Default is FALSE.
#' @param fname character string used to label file output. Default is "GWCoGAPS.AP.fixed"
#' @param PatternsMatchFN function to use for pattern matching across sets
#' @param Cut number of branches at which to cut dendrogram used in patternMatch4Parallel
#' @param minNS minimum of individual set contributions a cluster must contain
#' @param ... additional parameters to be fed into \code{gapsRun} and \code{gapsMapRun}
#' @export
#' @seealso \code{\link{gapsRun}}, \code{\link{patternMatch4Parallel}}, and \code{\link{gapsMapRun}}
#' @examples \dontrun{
#' GWCoGAPS(nCores=NA, D, S, nFactor, nSets,saveBySetResults=TRUE, fname=fname,
#' PatternsMatchFN = patternMatch4Parallel,numSnapshots=numSnapshots,minNS=minNS)
#' }
#'
GWCoGAPS <- function(D, S, params, nSets, nCores=NA, saveBySetResults=FALSE,
fname="GWCoGAPS.AP.fixed", PatternsMatchFN = patternMatch4Parallel, Cut=NA,
minNS=NA, ...)
{
    # v2 style parameters compatibility
    params <- oldParams(params, ...)

    # establish the number of cores that you are able to use
    if(is.na(nCores))
        nCores <- nSets
    registerDoParallel(cores=nCores)

    # break the data into sets
    genesInSets <- createGWCoGAPSSets(data=D, nSets=nSets, keep=FALSE)

    # set gene min to 0
    D <- sweep(D, 1, apply(D,1,function(x) pmin(0,min(x))))

    #generate seeds for parallelization
    nut <- generateSeeds(chains=nSets, seed=-1)

    # run CoGAPS for each set
    AP <- foreach(i=1:nSets) %dopar%
    {
        D <- as.matrix(D[genesInSets[[i]],])
        S <- as.matrix(S[genesInSets[[i]],])
        params <- setParam(params, 'seed', nut[i])
        CoGAPS(D=D, S=S, params)
    }

    BySet <- reOrderBySet(AP=AP, nFactor=nFactor, nSets=nSets)

    #run postpattern match function
    if (is.na(Cut))
        Cut <- nFactor
    matchedPs <- patternMatch4Parallel(Ptot=BySet$P, nP=nFactor, nSets=nSets,
        cnt=Cut, minNS=minNS, bySet=TRUE)

    #save BySet outputs
    class(AP) <- append(class(AP),"CoGAPS")
    if (saveBySetResults==TRUE)
    {
        save(AP,BySet,matchedPs,file=sprintf('APbySet.%s.nP%d.set%d.Rda',fname, nFactor,nSets))
        message(sprintf('APbySet.%s.nP%d.set%d.Rda',fname, nFactor,nSets))
    }

    PbySet <- matchedPs[["PBySet"]]
    matchedPs <- matchedPs[[1]]

    # generate seeds for parallelization
    nut <- generateSeeds(chains=nSets, seed=-1)

    # final number of factors
    nFactorFinal <- dim(matchedPs)[1]

    # run fixed CoGAPS
    params <- setParam('fixedMatrix', as.matrix(matchedPs))
    params <- setParam('whichMatrixFixed', 'P')
    params <- setParam('nFactor', nFactorFinal)
    Fixed <- foreach(i=1:nSets) %dopar%
    {
    	D <- as.matrix(D[genesInSets[[i]],])
    	S <- as.matrix(S[genesInSets[[i]],])
        params <- setParam(params, 'seed', nut[i])
        CoGAPS(D=D, S=S, params)
    }

    #extract A and Asds
    As4fixPs <- postFixed4Parallel(AP.fixed=Fixed, setPs=matchedPs)

    #save final
    AP.fixed <- list("Amean"=As4fixPs$A, "Asd"=As4fixPs$Asd, "Pmean"=matchedPs,
        "PbySet"=PbySet)
    class(AP.fixed) <- append(class(AP.fixed),"CoGAPS")
    save(AP.fixed, file=paste(fname, ".Rda", sep=""))
    message(paste(fname, ".Rda", sep=""))
}

#' Create Gene Sets for GWCoGAPS
#'
#'\code{createGWCoGAPSSets} factors whole genome data into randomly generated sets for indexing;
#'
#'@param data data matrix with unique rownames
#'@param nSets number of sets for parallelization
#'@param outRDA name of output file
#'@param keep logical indicating whether or not to save gene set list. Default is TRUE.
#'@export
#'@return list with randomly generated sets of genes from whole genome data
#'@examples \dontrun{
#'createGWCoGAPSSet(D,nSets=nSets)
#'}
createGWCoGAPSSets <- function(data, nSets, outRDA="GenesInCoGAPSSets.Rda",
keep=TRUE)
{
    genes <- rownames(data)
    setSize <- floor(length(genes) / nSets)
    genesInSets <- list()
    for (set in 1:nSets)
    {
        if (set != nSets)
            genesInSets[[set]] <- sample(genes,setSize)
        if (set == nSets)
            genesInSets[[set]] <- genes
        genes <- genes[!genes %in% genesInSets[[set]]]
    }

    if (!identical(sort(unlist(genesInSets)),sort(rownames(data))))
        warning("Gene identifiers not unique!")
    if (keep == TRUE)
        save(list=c('genesInSets'), file=outRDA)
    return(genesInSets)
}

#' Generate Seeds for Multiple Concurrent Runs
#'
#' @param chains number of seeds to generate (number of chains to run)
#' @param seed positive values are kept, negative values will be overwritten
#'  by a seed generated from the current time
#' @return vector of randomly generated seeds
#' @export
#' @examples
#' generateSeeds(chains=2, seed=-1)
generateSeeds <- function(chains=2, seed=-1)
{
    if (chains < 2 | (as.integer(chains) != chains))
        stop("chains must be >= 2 and an integer")

    if (seed < 0)
    {
        secs <- as.numeric(difftime(Sys.time(), paste(Sys.Date(), "00:00"),
            units="secs"))
        return(seq_len(chains) * round(secs))
    }
    else
    {
        return(seq_len(chains) * seed)
    }
}

#' reOrderBySet
#'
#' @description <restructures output of gapsRun into a list containing each sets solution for Amean, Pmean, and Asd>
#' @param AP output of gapsRun in parallel
#' @param nFactor number of patterns
#' @param nSets number of sets
#'
#' @return a list containing the \code{nSets} sets solution for Amean under "A", Pmean under "P", and Asd under "Asd"
#' @export
#'
#' @examples \dontrun{
#' reOrderBySet(AP,nFactor,nSets)
#' }
#'
reOrderBySet<-function(AP, nFactor, nSets)
{
    P <- do.call(rbind,lapply(AP, function(x) x$Pmean))
    rownames(P) <- paste(rep(1:nSets,each=nFactor),rep(1:nFactor,nSets),sep=".")
    A <- lapply(AP, function(x) x$Amean)
    Asd <- lapply(AP, function(x) x$Asd)
    names(A) <- names(Asd) <- paste(rep("Set", nSets), rep(1:nSets), sep="")
    return(list("A"=A, "Asd"=Asd, "P"=P))
}

#' postFixed4Parallel
#'
#' @param AP.fixed output of parallel gapsMapRun calls with same FP
#' @param setPs data.frame with rows giving fixed patterns for P used as input for gapsMapRun
#'
#' @return list of two data.frames containing the A matrix estimates or their corresponding standard deviations
#' from output of parallel gapsMapRun
#' @export
postFixed4Parallel <- function(AP.fixed=NA, setPs=NA)
{
    ASummary <- do.call(rbind,lapply(AP.fixed, function(x) x$Amean))
    Asd <- do.call(rbind,lapply(AP.fixed, function(x) x$Asd))
    #PSummary <- do.call(rbind,lapply(AP.fixed, function(x) x$Pmean))
    PSummary <- AP.fixed[[1]]$Pmean

    Pmax <- apply(PSummary,1,max)
    Pneu <- sweep(PSummary,1,Pmax,FUN="/")
    Aneu <- sweep(ASummary,2,Pmax,FUN="*")

    X <- apply(Pneu,1,range)
    Y <- apply(setPs,1,range)
    colnames(X) <- colnames(Y)
    if (!all.equal(X, Y, tolerance=1e-5))
        warning("Patterns do not match fixed values.")

    As4fixPs <- list("A"=Aneu, "Asd"=Asd)
    return(As4fixPs)
}
