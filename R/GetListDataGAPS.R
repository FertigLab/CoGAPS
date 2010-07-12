# formats raw table data into a format required for MCMC with JAGS
GetDataListGAPS <- function(data, unc, numPatterns,
                            MaxAtomsA, alphaA, lambdaA, 
                            MaxAtomsP, alphaP, lambdaP, SAIter, thin) {

  # convert to matrix for output
  D <- as.matrix(data)

  # process the uncertainty
  S <- as.matrix(unc)

  # check the dimensionality
  if (nrow(S) != nrow(D) || ncol(S) != ncol(D)) {
    stop(cat("Dimensions of data (",
             nrow(D), ", ", ncol(D),
             ") must match dimensions of uncertainty (",
             nrow(S), ", ", ncol(S), ").\n"))
  }

  # check that D and S are meet positivity requirements positive
  if (any(D < 0)) {
    warning("Data should be non-negative.\n");
  }

  if (any(S <= 0)) {
    stop("Uncertainty must be strictly positive.\n");
  }
    

  datalist = list(
    numGenes = nrow(data), numSamples = ncol(data),
    numPatterns = numPatterns,
    D = D, S = S,
    ADim = matrix(data=1,nrow=nrow(data), ncol=numPatterns),
    MaxAtomsA = MaxAtomsA, alphaA = alphaA, lambdaA = lambdaA,
    PDim = matrix(data=1,nrow=numPatterns, ncol=ncol(data)),
    MaxAtomsP = MaxAtomsP, alphaP = alphaP, lambdaP = lambdaP,
    SAIterA = SAIter, SAIterP = SAIter, fileID = thin)

  return(datalist)
  

  
}
