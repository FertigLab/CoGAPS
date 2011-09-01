# runs matrix decomposition with MCMC using rjags
RunGAPS <- function(data, unc, numPatterns,
                    MaxAtomsA, alphaA, lambdaA, 
                    MaxAtomsP, alphaP, lambdaP, 
                    SAIter, iter, thin=1) {


  # get the data needed for MCMC
  matrixDecompData <- GetDataListGAPS(data, unc, numPatterns,
                                      MaxAtomsA, alphaA, lambdaA,
                                      MaxAtomsP, alphaP, lambdaP,
                                      SAIter, thin)
  matrixDecompInits <- list(A = matrix(
                              data=0, matrixDecompData$numGenes, numPatterns),
                            P = matrix(
                              data=0, numPatterns, matrixDecompData$numSamples))

  # create the model file used by this decomposition
  FID <- sprintf('%.10f', thin - floor(thin))
  ModelFile <- paste("GAPS", FID, "bug", sep=".")

  cat("var\n", file=ModelFile)
  cat("\tA[numGenes, numPatterns], ADim[numGenes, numPatterns],\n",
      file=ModelFile, append=TRUE)
  cat("\tP[numPatterns, numSamples], PDim[numPatterns, numSamples],\n",
      file=ModelFile, append=TRUE)
  cat("\tD[numGenes, numSamples], S[numGenes, numSamples];\n",
      file=ModelFile, append=TRUE)
  cat("\nmodel {\n", file=ModelFile, append=TRUE)
  cat("A[,] ~ datomicgibbs(ADim, MaxAtomsA, alphaA, lambdaA, SAIterA, fileID)\n", file=ModelFile, append=TRUE)
  cat("P[,] ~ datomicgibbs(PDim, MaxAtomsP, alphaP, lambdaP, SAIterP, fileID)\n", file=ModelFile, append=TRUE)
  cat("D[,] ~ gapsnorm(A,P,S)\n", file=ModelFile, append=TRUE)
  cat("}\n", file=ModelFile, append=TRUE) 

  # set the model being used
  matrixDecompModel <- jags.model(file=ModelFile,
                                  data=matrixDecompData,
                                  inits=matrixDecompInits,
                                  n.adapt=SAIter)
  # update and specify which nodes are monitored
  matrixDecompOutput <- jags.samples(model=matrixDecompModel,
                                     variable.names=c('A','P'),  
                                     n.iter=iter,
                                     thin=iter)

  file.remove(ModelFile)
  
  return(matrixDecompOutput)
  

}

