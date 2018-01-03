
# gapsRun: function to call C++ cogaps code
# History: v1.0 CK with MFO edits, August 2014

# Inputs: D - data matrix
#         S - uncertainty matrix (std devs for chi-squared of Log Likelihood)
#         nFactor - number of patterns (basis vectors, metagenes)
#         simulation_id - name to attach to atoms files if created
#         nEquil - number of iterations for burn-in
#         nSample - number of iterations for sampling
#         nOutR - how often to print status into R by iterations
#         output_atomic - whether to write atom files (large)
#         alphaA, alphaP - sparsity parameters for A and P domains
#         max_gibbmass_paraA(P) - limit truncated normal to max size
#         nMaxA, nMaxP - PRESENTLY UNUSED, future = limit number of atoms


# Output: list with A and P matrix estimates, chi-squared and atom
#         numbers of sample by iteration, and chi-squared of mean
#'\code{gapsRun} calls the C++ MCMC code and performs Bayesian
#'matrix factorization returning the two matrices that reconstruct
#'the data matrix
#'
#'@param D data matrix
#'@param S uncertainty matrix (std devs for chi-squared of Log Likelihood)
#'@param ABins a matrix of same size as A which gives relative
#' probability of that element being non-zero
#'@param PBins a matrix of same size as P which gives relative
#' probability of that element being non-zero
#'@param nFactor number of patterns (basis vectors, metagenes), which must be
#'greater than or equal to the number of rows of FP
#'@param  simulation_id name to attach to atoms files if created
#'@param  nEquil number of iterations for burn-in
#'@param  nSample number of iterations for sampling
#'@param  nOutR how often to print status into R by iterations
#'@param  output_atomic whether to write atom files (large)
#'@param fixedBinProbs Boolean for using relative probabilities
#' given in Abins and Pbins
#'@param fixedDomain character to indicate whether A or P is
#' domain for relative probabilities
#'@param numSnapshots the number of individual samples to capture
#'@param alphaA sparsity parameter for A domain
#'@param nMaxA PRESENTLY UNUSED, future = limit number of atoms
#'@param max_gibbmass_paraA limit truncated normal to max size
#'@param alphaP sparsity parameter for P domain
#'@param nMaxP PRESENTLY UNUSED, future = limit number of atoms
#'@param max_gibbmass_paraP limit truncated normal to max size
#'@param seed Set seed for reproducibility. Positive values provide initial seed, negative values just use the time.
#'@param messages Display progress messages
#'@param singleCellRNASeq T/F indicating if the data is single cell RNA-seq data
#'@param fixedPatterns matrix of fixed values in either A or P matrix
#'@param whichMatrixFixed character to indicate whether A or P matrix
#'  contains the fixed patterns
#'@export

#--CHANGES 1/20/15--
#Added FixedPatt frame to C++ version to match Ondrej's code for updating bin sizes based on an input matrix
gapsRun <- function(D, S, ABins = data.frame(), PBins = data.frame(),
                    nFactor = 7, simulation_id = "simulation",
                    nEquil = 1000, nSample = 1000, nOutR = 1000,
                    output_atomic = FALSE, fixedBinProbs = FALSE,
                    fixedDomain = "N", numSnapshots = 0, alphaA = 0.01,
                    nMaxA = 100000, max_gibbmass_paraA = 100.0,
                    alphaP = 0.01, nMaxP = 100000, max_gibbmass_paraP = 100.0,
                    seed=-1, messages=TRUE, singleCellRNASeq=FALSE,
                    fixedPatterns = matrix(0), whichMatrixFixed = 'N')
{
    # Floor the parameters that are integers to prevent allowing doubles.
    nFactor <- floor(nFactor)
    nEquil <- floor(nEquil)
    nSample <- floor(nSample)
    nOutR <- floor(nOutR)
    numSnapshots <- floor(numSnapshots)
    nMaxA <- floor(nMaxA)
    nMaxP <- floor(nMaxP)

    # check the fixed patterns
    if ((whichMatrixFixed == 'A' & nrow(fixedPatterns) != nrow(D))
    | (whichMatrixFixed == 'P' & nrow(fixedPatterns) > nFactor)) 
    {
        stop('invalid number of rows in fixed pattern matrix')
    }
    if ((whichMatrixFixed == 'A' & ncol(fixedPatterns) > nFactor)
    | (whichMatrixFixed == 'P' & ncol(fixedPatterns) != ncol(D))) 
    {
        stop('invalid number of cols in fixed pattern matrix')
    }

    geneNames <- rownames(D);
    sampleNames <- colnames(D);

    # label patterns as Patt N
    patternNames <- c("0");
    for(i in 1:nFactor)
    {
        patternNames[i] <- paste('Patt', i);
    }

    # call to C++ Rcpp code
    cogapResult <- cogaps(as.matrix(D), as.matrix(S), nFactor, alphaA, alphaP,
        nEquil, floor(nEquil/10), nSample, max_gibbmass_paraA, max_gibbmass_paraP,
        fixedPatterns, whichMatrixFixed, seed, messages, singleCellRNASeq,
        nOutR, numSnapshots)

    # convert returned files to matrices to simplify visualization and processing
    cogapResult$Amean <- as.matrix(cogapResult$Amean);
    cogapResult$Asd <- as.matrix(cogapResult$Asd);
    cogapResult$Pmean <- as.matrix(cogapResult$Pmean);
    cogapResult$Psd <- as.matrix(cogapResult$Psd);

    ## label matrices
    colnames(cogapResult$Amean) <- patternNames;
    rownames(cogapResult$Amean) <- geneNames;
    colnames(cogapResult$Asd) <- patternNames;
    rownames(cogapResult$Asd) <- geneNames;
    colnames(cogapResult$Pmean) <- sampleNames;
    rownames(cogapResult$Pmean) <- patternNames;
    colnames(cogapResult$Psd) <- sampleNames;
    rownames(cogapResult$Psd) <- patternNames;

    ## calculate chi-squared of mean, this should be smaller than individual
    ## chi-squared sample values if sampling is good
    calcChiSq <- c(0);
    MMatrix <- (cogapResult$Amean %*% cogapResult$Pmean);

    for(i in 1:(nrow(MMatrix)))
    {
        for(j in 1:(ncol(MMatrix)))
        {
            calcChiSq <- calcChiSq + ((D[i,j] - MMatrix[i,j])/S[i,j])^2;
        }
    }

    cogapResult = c(cogapResult, calcChiSq);
    
    # names(cogapResult)[13] <- "meanChi2";

    message(paste("Chi-Squared of Mean:", calcChiSq))
    return(cogapResult);
}
