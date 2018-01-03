
# gapsMapRun: function to call C++ cogaps code
# History: v1.0 CK with MFO edits, August 2014

# Inputs: D - data matrix
#         S - uncertainty matrix (std devs for chi-squared of Log Likelihood)
#         FP - fixed A columns or P rows
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

#'\code{gapsMapRun} calls the C++ MCMC code and performs Bayesian
#'matrix factorization returning the two matrices that reconstruct
#'the data matrix; as opposed to gapsRun, this method takes an
#'additional input specifying set patterns in the P matrix
#'
#'@param D data matrix
#'@param S uncertainty matrix (std devs for chi-squared of Log Likelihood)
#'@param FP data.frame with rows giving fixed patterns for P
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
#'@param  fixedMatrix character indicating whether A or P matrix
#' has fixed columns or rows respectively
#'@param fixedBinProbs Boolean for using relative probabilities
#' given in Abins and Pbins
#'@param fixedDomain character to indicate whether A or P is
#' domain for relative probabilities
#'@param sampleSnapshots Boolean to indicate whether to capture
#' individual samples from Markov chain during sampling
#'@param numSnapshots the number of individual samples to capture
#'@param alphaA sparsity parameter for A domain
#'@param nMaxA PRESENTLY UNUSED, future = limit number of atoms
#'@param max_gibbmass_paraA limit truncated normal to max size
#'@param alphaP sparsity parameter for P domain
#'@param nMaxP PRESENTLY UNUSED, future = limit number of atoms
#'@param max_gibbmass_paraP limit truncated normal to max size
#'@param seed Set seed for reproducibility. Positive values provide initial seed, negative values just use the time.
#'@param messages Display progress messages
#'@export

gapsMapRun <- function(D, S, FP, ABins = data.frame(), PBins = data.frame(), nFactor = 7, simulation_id = "simulation",
                       nEquil = 1000, nSample = 1000, nOutR = 1000, output_atomic = FALSE, fixedMatrix = "P",
                       fixedBinProbs = FALSE, fixedDomain = "N", sampleSnapshots = TRUE, numSnapshots = 100, alphaA = 0.01,
                       nMaxA = 100000, max_gibbmass_paraA = 100.0, alphaP = 0.01, nMaxP = 100000, max_gibbmass_paraP = 100.0,
                       seed=-1, messages=TRUE)
{
    gapsRun(D, S, ABins, PBins, nFactor, simulation_id, nEquil, nSample, nOutR, output_atomic, fixedBinProbs,
        fixedDomain, sampleSnapshots, numSnapshots, alphaA, nMaxA, max_gibbmass_paraA, alphaP, nMaxP,
        max_gibbmass_paraP, seed, messages, FALSE, FP, 'P')
}
