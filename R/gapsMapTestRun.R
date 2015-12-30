
# gapsMapTestRun: function to call C++ cogaps code
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

#'\code{gapsMapTestRun} calls the C++ MCMC code and performs Bayesian
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
#'@param alphaA sparsity parameter for A domain
#'@param nMaxA PRESENTLY UNUSED, future = limit number of atoms
#'@param max_gibbmass_paraA limit truncated normal to max size
#'@param alphaP sparsity parameter for P domain
#'@param nMaxP PRESENTLY UNUSED, future = limit number of atoms
#'@param max_gibbmass_paraP limit truncated normal to max size
#'@export

gapsMapTestRun <- function(D, S, FP, ABins = data.frame(), PBins = data.frame(),
                           nFactor = 7, simulation_id = "simulation",
                           nEquil = 1000, nSample = 1000, nOutR = 1000,
                           output_atomic = FALSE, fixedMatrix = "P",
                           fixedBinProbs = FALSE, fixedDomain = "N",
                           alphaA = 0.01,  nMaxA = 100000,
                           max_gibbmass_paraA = 100.0, alphaP = 0.01,
                           nMaxP = 100000, max_gibbmass_paraP = 100.0)
{

  #Begin data type error checking code
  charDataErrors = c(!is.character(simulation_id), !is.character(fixedDomain), !is.character(fixedMatrix))
  charCheck = c("simulation_id", "fixedDomain", "fixedMatrix")

  boolDataErrors = c(!is.logical(output_atomic), !is.logical(fixedBinProbs))
  boolCheck = c("output_atomic", "fixedBinProbs")

  numericDataErrors = c(!is.numeric(nFactor), !is.numeric(nEquil), !is.numeric(nSample), !is.numeric(nOutR),
                        !is.numeric(alphaA), !is.numeric(nMaxA), !is.numeric(max_gibbmass_paraA), !is.numeric(alphaP),
                        !is.numeric(nMaxP), !is.numeric(max_gibbmass_paraP))
  numericCheck = c("nFactor", "nEquil", "nSample", "nOutR", "alphaA", "nMaxA",
                   "max_gibbmass_paraA", "alphaP",    "nMaxP", "max_gibbmass_paraP")

  dataFrameErrors = c(!is.data.frame(D), !is.data.frame(S), !is.data.frame(ABins), !is.data.frame(PBins), !is.data.frame(FP))
  dataFrameCheck = c("D", "S", "ABins", "PBins", "FP")

  matrixErrors = c(!is.matrix(D), !is.matrix(S), !is.matrix(ABins), !is.matrix(PBins), !is.matrix(FP))

  if(any(charDataErrors))
  {
    #Check which of these is not a string
    for(i in 1:length(charDataErrors))
    {
      if(charDataErrors[i] == TRUE)
      {
        stop(paste("Error in gapsRun: Argument",charCheck[i],"is of the incorrect type. Please see documentation for details."))
      }
    }
  }

  if(any(boolDataErrors))
  {
    #Check which of these is not a boolean
    for(i in 1:length(boolDataErrors))
    {
      if(boolDataErrors[i] == TRUE)
      {
        stop(paste("Error in gapsRun: Argument",boolCheck[i],"is of the incorrect type. Please see documentation for details."))
      }
    }
  }

  if(any(numericDataErrors))
  {
    #Check which of these is not an integer/double
    for(i in 1:length(numericDataErrors))
    {
      if(numericDataErrors[i] == TRUE)
      {
        stop(paste("Error in gapsRun: Argument",numericCheck[i],"is of the incorrect type. Please see documentation for details."))

      }
    }
  }


  #At least one of A, P, ABins, PBins, or FP is not a matrix or data.frame
  if(any((dataFrameErrors[1] && matrixErrors[1]), (dataFrameErrors[2] && matrixErrors[2]), (dataFrameErrors[3] && matrixErrors[3]), (dataFrameErrors[4] && matrixErrors[4]), (dataFrameErrors[5] && matrixErrors[5])))
  {
    for(i in 1:length(dataFrameCheck))
    {
      if((dataFrameErrors[i] && matrixErrors[i]) == TRUE)
      {
        stop(paste("Error in gapsRun: Argument", dataFrameCheck[i], "is not a matrix or data.frame. Please see documentation for details."))

      }
    }
  }

  #Floor the parameters that are integers to prevent allowing doubles.
  nFactor = floor(nFactor)
  nEquil = floor(nEquil)
  nSample = floor(nSample)
  nOutR = floor(nOutR)
  nMaxA = floor(nMaxA)
  nMaxP = floor(nMaxP)

  # pass all settings to C++ within a list
  #    if (is.null(P)) {
  Config = c(simulation_id, output_atomic, fixedBinProbs, fixedDomain, fixedMatrix);

  ConfigNums = c(nFactor, nEquil, nSample, nOutR, alphaA, nMaxA, max_gibbmass_paraA,
                 alphaP, nMaxP, max_gibbmass_paraP);

  #Begin logic error checking code

  #Check for negative or zero arguments
  if(any(ConfigNums <= 0))
  {
    stop("Error in gapsRun: Numeric Arguments cannot be non-zero!")
  }

  if((nOutR > nEquil) || (nOutR > nSample))
  {
    stop("Error in gapsRun: Cannot have more output steps than equilibration and/or sampling iterations.")
  }

  if(ncol(FP) != ncol(D))
  {
    stop("Error in gapsRun: Columns of Data Matrix and Fixed Pattern Matrix do not line up. Please see documentation for details.")
  }

  if(nFactor < (nrow(FP)))
  {
    stop("Error in gapsRun: Number of patterns cannot be less than the rows of the patterns to fix (FP). Please see documentation for details.")
  }

  if(nFactor > (ncol(D)))
  {
    warning("Warning in gapsRun: Number of requested patterns greater than columns of Data Matrix.")
  }

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
  cogapResult = cogapsMapTest(D, S, FP, ABins, PBins, Config, ConfigNums);

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
  names(cogapResult)[20] = "meanChi2";

  message(paste("Chi-Squared of Mean:",calcChiSq))

  return(cogapResult);
}
