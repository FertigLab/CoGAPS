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
#         lambdaA(P)_scale_factor - lambda factor in penalized likelihood

# Output: list with matrices, mean chi-squared value, and gene set
#         results

#'\code{CoGAPS} calls the C++ MCMC code through gapsRun and performs Bayesian
#'matrix factorization returning the two matrices that reconstruct
#'the data matrix and then calls calcCoGAPSStat to estimate gene set
#'activity with nPerm set to 500
#'
#'@param data data matrix
#'@param unc uncertainty matrix (std devs for chi-squared of Log Likelihood)
#'@param GStoGenes data.frame or list with gene sets
#'@param nFactor number of patterns (basis vectors, metagenes)
#'@param simulation_id name to attach to atoms files if created
#'@param plot logical to determine if plots produced
#'@param nPerm number of permutations for gene set test
#'@param nEquil number of iterations for burn-in
#'@param nSample number of iterations for sampling
#'@param nOutR how often to print status into R by iterations
#'@param output_atomic whether to write atom files (large)
#'@param alphaA sparsity parameter for A domain 
#'@param alphaP sparsity parameter for P domain 
#'@param max_gibbmass_paraA limit truncated normal to max size for A
#'@param max_gibbmass_paraP limit truncated normal to max size for P
#'@param nMaxA PRESENTLY UNUSED, future = limit number of atoms for A
#'@param nMaxP PRESENTLY UNUSED, future = limit number of atoms for P
#'@param lambdaA_scale_factor lambda factor in penalized likelihood for A
#'@param lambdaP_scale_factor lambda factor in penalized likelihood for P
#'@export

CoGAPS <- function(data, unc, GStoGenes, nFactor = "7", nEquil=1000,
                nSample=1000, nOutR=1000, output_atomic="false",
                simulation_id="simulation", plot=TRUE,nPerm=500,
                 alphaA = "0.01",  nMaxA = "100000", 
                 max_gibbmass_paraA = "100.0", lambdaA_scale_factor = "1.0", 
                 alphaP = "0.01", nMaxP = "100000", 
                 max_gibbmass_paraP = "100.0", lambdaP_scale_factor = "1.0") {


  # decompose the data
  matrixDecomp <- gapsRun(data, unc, nFactor=nFactor, 
    simulation_id=simulation_id, nEquil=nEquil, 
    nSample=nSample, nOutR = nOutR, output_atomic = output_atomic, 
    alphaA = alphaA,  nMaxA = nMaxA, 
    max_gibbmass_paraA = max_gibbmass_paraA, 
    lambdaA_scale_factor = lambdaA_scale_factor, 
    alphaP = alphaP, nMaxP = nMaxP, 
    max_gibbmass_paraP = max_gibbmass_paraP, 
    lambdaP_scale_factor = lambdaP_scale_factor)

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
