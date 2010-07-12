# run the full algorithm
CoGAPS <- function(data, unc, GStoGenes,
                 outputDir, outputBase="",
                 sep = "\t", isPercentError=FALSE,
                 numPatterns,
                 MaxAtomsA=2^32, alphaA=0.01,
                 MaxAtomsP=2^32, alphaP=0.01,
                 SAIter=1000000000, iter = 500000000, thin=-1,
                 nPerm=500, verbose=TRUE, plot=FALSE, keepChain=FALSE) {


  # decompose the data
  matrixDecomp <- GAPS(data=data, unc=unc,
                       outputDir=outputDir, outputBase=outputBase, 
		       sep=sep, isPercentError=isPercentError,
                       numPatterns=numPatterns, MaxAtomsA=MaxAtomsA, 
		       alphaA=alphaA, MaxAtomsP=MaxAtomsP,
                       alphaP=alphaP, SAIter=SAIter, iter=iter, thin=thin,
                       verbose=verbose, keepChain=keepChain)

  # plot patterns and show heatmap of Anorm
  if (plot) {
    plotGAPS(matrixDecomp$Amean, matrixDecomp$Pmean)
  }

  # compute the gene set scores
  GSP <- calcCoGAPSStat(matrixDecomp$Amean, matrixDecomp$Asd, GStoGenes, nPerm)

  # output gene set results
  if (outputDir != "") {
     outputGSPBase <- paste(outputDir, outputBase, sep="/")
  } else {
     outputGSPBase <- outputBase
  }
  if (outputBase != "") {
    outputGSPBase <- paste(outputGSPBase, "_", sep="")
  }
  
  UpFile <- paste(outputGSPBase, "GSUpStat.txt", sep="")
  DownFile <- paste(outputGSPBase, "GSDownStat.txt", sep="")
  ActEstFile <- paste(outputGSPBase, "GSActEst.txt", sep="")

  write.table(GSP$GSUpreg, file=UpFile, sep=sep)
  write.table(GSP$GSDownreg, file=DownFile, sep=sep)
  write.table(GSP$GSActEst, file=ActEstFile, sep=sep)

  return(list(meanChi2=matrixDecomp$meanChi2,
              D=matrixDecomp$D, Sigma=matrixDecomp$Sigma,
              Amean=matrixDecomp$Amean, Asd=matrixDecomp$Asd,
              Pmean=matrixDecomp$Pmean, Psd=matrixDecomp$Psd,
              meanMock=matrixDecomp$meanMock,
              GSUpreg=GSP$GSUpreg, GSDownreg=GSP$GSDownreg, GSActEst=GSP$GSActEst))
  
  
}
