
# gapsRun: function to call C++ cogaps code
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
#         lambdaA(P)_scale_factor - lambda factor in penalized likelihood

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
#'@param nFactor number of patterns (basis vectors, metagenes), which must be
#'greater than or equal to the number of rows of FP
#'@param  simulation_id name to attach to atoms files if created
#'@param  nEquil number of iterations for burn-in
#'@param  nSample number of iterations for sampling
#'@param  nOutR how often to print status into R by iterations
#'@param  output_atomic whether to write atom files (large)
#'@param  alphaA sparsity parameter for A domain 
#'@param  alphaP sparsity parameter for P domain 
#'@param  max_gibbmass_paraA limit truncated normal to max size in A
#'@param  max_gibbmass_paraP limit truncated normal to max size in P
#'@param  nMaxA PRESENTLY UNUSED, future = limit number of atoms in A
#'@param  nMaxP PRESENTLY UNUSED, future = limit number of atoms in P
#'@param  lambdaA_scale_factor lambda factor in penalized likelihood in A
#'@param  lambdaP_scale_factor lambda factor in penalized likelihood in P
#'@export

gapsMapRun <- function(D, S, FP, nFactor = "7", simulation_id = "simulation", nEquil = "1000", nSample = "1000", nOutR = 1000, output_atomic = "FALSE", alphaA = "0.01",  nMaxA = "100000", max_gibbmass_paraA = "100.0", lambdaA_scale_factor = "1.0", alphaP = "0.01", nMaxP = "100000", max_gibbmass_paraP = "100.0", lambdaP_scale_factor = "1.0")
{
    
    # scale fixed patterns to have norm 1
    FPS <- apply(FP,1,sum)
    if (any(FPS==0) | any(FP<0)) {
        stop('Fixed patterns in FP must be non-negative and have non-zero norm')
    }
    FP <- sweep(FP,1,FPS,FUN="/")
    
    # pass all settings to C++ within a list
    #    if (is.null(P)) {
        Config = c(nFactor, simulation_id, nEquil, nSample, nOutR,
        output_atomic, alphaA, nMaxA, max_gibbmass_paraA, lambdaA_scale_factor,
        alphaP, nMaxP, max_gibbmass_paraP, lambdaP_scale_factor)
        #        P <- as.data.frame(matrix(nrow=1,c(1,1,1))) # make something to pass
        #    } else {
        #        Config = c(nFactor, simulation_id, nEquil, nSample, nOutR,
        #        output_atomic, alphaA, nMaxA, max_gibbmass_paraA, lambdaA_scale_factor,
        #        alphaP, nMaxP, max_gibbmass_paraP, lambdaP_scale_factor, 1)
 
 #   }
	
	geneNames = rownames(D);
	sampleNames = colnames(D);
	
    # label patterns as Patt N
	patternNames = c("0");
	for(i in 1:nFactor)
	{
		patternNames[i] = paste('Patt', i);
	}

    # call to C++ Rcpp code
	cogapResult = cogapsmap(D, S, FP, Config);

    # convert returned files to matrices to simplify visualization and processing
	cogapResult$Amean = as.matrix(cogapResult$Amean);
	cogapResult$Asd = as.matrix(cogapResult$Asd);
	cogapResult$Pmean = as.matrix(cogapResult$Pmean);
	cogapResult$Psd = as.matrix(cogapResult$Psd);

    # label matrices
	colnames(cogapResult$Amean) = patternNames;
	rownames(cogapResult$Amean) = geneNames;
	colnames(cogapResult$Asd) = patternNames;
	rownames(cogapResult$Asd) = geneNames;
	colnames(cogapResult$Pmean) = sampleNames;
	rownames(cogapResult$Pmean) = patternNames;
	colnames(cogapResult$Psd) = sampleNames;
	rownames(cogapResult$Psd) = patternNames;

    # calculate chi-squared of mean, this should be smaller than individual
    # chi-squared sample values if sampling is good
	calcChiSq = c(0);
	MMatrix = (cogapResult$Amean %*% cogapResult$Pmean);

	
	for(i in 1:(nrow(MMatrix)))
	{
		for(j in 1:(ncol(MMatrix)))
		{
			calcChiSq = (calcChiSq) + ((D[i,j] - MMatrix[i,j])/(S[i,j]))^(2);
		}
	}
	
	cogapResult = c(cogapResult, calcChiSq);
	names(cogapResult)[10] = "meanChi2";
    
    print(paste("Chi-Squared of Mean:",calcChiSq))
	
	return(cogapResult);
}
