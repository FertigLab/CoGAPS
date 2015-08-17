# CoGAPS: function to call gapsRun and calculate gene set
#         statistics on result, also to produce A and P plots
# History: v1.0 - EJF original CoGAPS, MFO modified to to match
#                 new variable names

# Inputs: data - data matrix
#         unc - uncertainty matrix (std devs for chi-squared of Log Likelihood)
#         GStoGenes - data.frame or list of gene sets
#         nFactor - number of patterns (basis vectors, metagenes) 
#         nEquil - number of iterations for burn-in
#         nSample - number of iterations for sampling
#         nOutR - how often to print status into R by iterations
#         output_atomic - whether to write atom files (large)
#         simulation_id - name to attach to atoms files if created
#         plot - logical to determine if plots produced
#         nPerm - number of permutations for gene set test
#         alphaA, alphaP - sparsity parameters for A and P domains
#         max_gibbmass_paraA(P) - limit truncated normal to max size
#         nMaxA, nMaxP - PRESENTLY UNUSED, future = limit number of atoms


# Output: list with matrices, mean chi-squared value, and gene set
#         results

#'\code{CoGAPS} calls the C++ MCMC code through gapsRun and performs Bayesian
#'matrix factorization returning the two matrices that reconstruct
#'the data matrix and then calls calcCoGAPSStat to estimate gene set
#'activity with nPerm set to 500
#'
#'@param data data matrix
#'@param unc uncertainty matrix (std devs for chi-squared of Log Likelihood)
#'@param ABins a matrix of same size as A which gives relative
#' probability of that element being non-zero
#'@param PBins a matrix of same size as P which gives relative
#' probability of that element being non-zero
#'@param GStoGenes data.frame or list with gene sets
#'@param nFactor number of patterns (basis vectors, metagenes)
#'@param simulation_id name to attach to atoms files if created
#'@param nEquil number of iterations for burn-in
#'@param nSample number of iterations for sampling
#'@param nOutR how often to print status into R by iterations
#'@param output_atomic whether to write atom files (large)
#'@param fixedBinProbs Boolean for using relative probabilities
#' given in Abins and Pbins
#'@param fixedDomain character to indicate whether A or P is
#' domain for relative probabilities
#'@param sampleSnapshots Boolean to indicate whether to capture
#' individual samples from Markov chain during sampling
#'@param numSnapshots the number of individual samples to capture
#'@param plot Boolean to indicate whether to produce output graphics
#'@param nPerm number of permutations in gene set test
#'@param alphaA sparsity parameter for A domain
#'@param nMaxA PRESENTLY UNUSED, future = limit number of atoms
#'@param alphaP sparsity parameter for P domain
#'@param max_gibbmass_paraA limit truncated normal to max size
#'@param nMaxP PRESENTLY UNUSED, future = limit number of atoms
#'@param max_gibbmass_paraP limit truncated normal to max size
#'@export

CoGAPS <- function(data, unc, ABins = data.frame(), PBins = data.frame(), GStoGenes, nFactor = 7, simulation_id="simulation", nEquil=1000,
                nSample=1000, nOutR=1000, output_atomic=FALSE,
                fixedBinProbs = FALSE, fixedDomain = "N",
                sampleSnapshots = TRUE, numSnapshots = 100,
                plot=TRUE, nPerm=500,
                alphaA = 0.01,  nMaxA = 100000,
                max_gibbmass_paraA = 100.0,
                alphaP = 0.01, nMaxP = 100000,
                max_gibbmass_paraP = 100.0) {


  # decompose the data
  matrixDecomp <- gapsRun(data, unc, ABins, PBins, nFactor, simulation_id,
					nEquil, nSample, nOutR, output_atomic, fixedBinProbs, 
					fixedDomain, sampleSnapshots, numSnapshots, alphaA,  nMaxA, max_gibbmass_paraA, 
					alphaP, nMaxP, max_gibbmass_paraP)

  # plot patterns and show heatmap of Anorm
  if (plot) {
    plotGAPS(matrixDecomp$Amean, matrixDecomp$Pmean)
  }

  # compute the gene set scores
  GSP <- calcCoGAPSStat(matrixDecomp$Amean, matrixDecomp$Asd, GStoGenes, nPerm)

  return(list(meanChi2=matrixDecomp$calcChiSq,
              D=data, Sigma=unc,
              Amean=matrixDecomp$Amean, Asd=matrixDecomp$Asd,
              Pmean=matrixDecomp$Pmean, Psd=matrixDecomp$Psd,
              GSUpreg=GSP$GSUpreg, GSDownreg=GSP$GSDownreg, GSActEst=GSP$GSActEst))
  
  
}
