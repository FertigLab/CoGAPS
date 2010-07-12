# driver for the matrix decomposition
GAPS <- function(data, unc,
                 outputDir, outputBase="",
                 sep = "\t", isPercentError=FALSE,
                 numPatterns,
                 MaxAtomsA=2^32, alphaA=0.01,
                 MaxAtomsP=2^32, alphaP=0.01,
                 SAIter=1000000000, iter = 500000000, thin=-1,
                 verbose=TRUE, keepChain = FALSE) {
  
  # formatting data for use by algorithm
  if (verbose) {
    message("Formatting data for GAPS");
  }

  # process data
  dataInFile = tolower(class(data))=="character"

  # read in data from a file
  if (dataInFile) {

    if (verbose) {
      message(paste("Reading in data from", data, sep=" "))
    }

    D <- as.matrix(read.table(file = data, sep = sep,
                              row.names=1, header=TRUE))
    
  } else {
    # data is input in matrix / array form
    D <- as.matrix(data)
  }

  if (numPatterns >= nrow(D) || numPatterns >= ncol(D)) {
    stop(sprintf("Cannot decompose data into more patterns (%d) than rows (%d) or columns (%d).", numPatterns, nrow(D), ncol(D)))
  }
  
  # process uncertainty
  sdInFile = tolower(class(unc))=="character"
  if (sdInFile) {

    # uncertainty is specified in a file
    if (isPercentError) {
      stop(paste("Uncertainty specified in file", unc, 
                  "must be specified as additive error."))
    }

    if (verbose) {
       message(paste("Reading in uncertainty from", unc, sep=" "))
    }
    
    Sigma <- as.matrix(read.table(file = unc, sep = sep,
                                  row.names=1,header=TRUE))
    
  } else {

    
    nrunc <- nrow(unc)
    ncunc <- ncol(unc)
    if (length(unc)==1) {
      nrunc <- 1
      ncunc <- 1
    }
    
    # see if input uncertainty is a scalar
    if ( nrunc == 1 && ncunc == 1 ) {

      if (verbose) {
        message("Using single unc value specified in input.")
      }

      # compute the error as either fractional or additive
      if (isPercentError) {
        Sigma <- pmax(unc*abs(D),unc)
      } else {
        Sigma <- matrix(data=unc, nrow=nrow(D), ncol=ncol(D))
      }

      # check that the uncertainty is positive
      if (any(Sigma <= 0)) {
        stop(paste("Uncertainty must be positive (Sigma = ", Sigma, ").", sep=""))
      }
    } else {

      if (verbose) {
        message("Using matrix unc value specified in input.")
      }

      # check that dimensions of uncertainty and data match
      if (nrunc != nrow(D) || ncunc != ncol(D)) {
        stop(paste("Dimensions of data (", nrow(D), ", ", ncol(D),
                    ") and uncertainty (", nrow(Sigma), ", ", ncol(Sigma), ") do not agree",
                    sep = ""))
      }
      # check that the uncertainty is positive
      if (any(unc <= 0)) {
        stop(paste("Uncertainty must be positive (unc = ", unc, ").", sep=""))
      }

      # ensure that matrix errors additive                                                        
      if (isPercentError) {
        stop("Uncertainty specified in an array must be additive error.")
      }
      Sigma <- as.matrix(unc)

    }

  }

  # make a temporary file to outputBase which contains the row and column names
  if (outputDir != "" && outputDir != ".") {
     fullOutputPath <- paste(outputDir,outputBase,sep="/")
     if (!file.exists(outputDir)) {
        dir.create(outputDir)
     }
  } else {
     fullOutputPath <- outputBase
  }
  

  if (verbose) {
    message(paste("Decomposing data of size:",
                nrow(D), "by", ncol(D),
                "into", numPatterns, "patterns.", sep=" "))
  } 
  
  # get expected mass of each matrix element based on data
  if (verbose) {
    message("Initializing distribution parameters")
  }

  # computing the mean mass of atoms in A and P based on the data matrix
  lambdaA <- alphaA * sqrt(numPatterns) / sqrt(mean(D)) 
  lambdaP <- alphaP * sqrt(numPatterns) / sqrt(mean(D))

  # estimate number of outputs based on expected number of matrix elements
  if (thin <= 0) {
    if (iter < 10000) {
      outThin <- 1
    } else {
      outThin <- iter/10000
    }
  } else {
    outThin <- thin
  }

  # assign the unique file identifier to be the number after the decimal
  if (floor(outThin) - outThin < 1.e-10) {
    outThin <- outThin + runif(1)
  }


  # initialize output arrays, including chain dimensions
  Amean   <- array(0., c(nrow(D), numPatterns))
  row.names(Amean) <- row.names(D)
  colnames(Amean) <- paste("p",1:numPatterns,sep="")
  Asd     <- array(0., c(nrow(D), numPatterns))
  row.names(Asd) <- row.names(D)
  colnames(Amean) <- paste("p",1:numPatterns,sep="")
  Pmean   <- array(0., c(numPatterns, ncol(D)))
  row.names(Pmean) <-   paste("p",1:numPatterns,sep="")
  colnames(Pmean) <- colnames(D)
  Psd     <- array(0., c(numPatterns, ncol(D)))
  row.names(Psd) <- paste("p",1:numPatterns,sep="")
  colnames(Psd) <- colnames(D)

  meanChi2 <- 0. 
  meanMock <- array(0., c(nrow(D),ncol(D)))

  
  if (verbose) {
    message(paste("Running GAPS."))
    message(paste(SAIter, "burn-in steps;", iter, "steps.", sep=" "))
  }
    
  # run the mcmc
  mdout <- RunGAPS(D, Sigma, numPatterns,
                   MaxAtomsA, alphaA, lambdaA,
                   MaxAtomsP, alphaP, lambdaP,
                   SAIter, iter, outThin)

  if (verbose) {
    message(paste("Normalizing results."))
  }

  # obtain the mean and standard deviation of A and P from the chain files
  FID <- sprintf('%.10f', outThin - floor(outThin))
  AFile <- paste("AResults", FID, ".ChainValues.txt", sep="")
  PFile <- paste("PResults", FID, ".ChainValues.txt", sep="")
  normmd <- AtomChain2APNorm(Afile      = AFile,
                             Pfile      = PFile,
                             nGene      = nrow(D),
                             nPattern   = numPatterns,
                             nSample    = ncol(D))
    
  # extract chain data and store in above arrays
  Amean[,] <- normmd$Amean
  Asd[,] <- normmd$Asd
  Pmean[,] <- normmd$Pmean
  Psd[,] <- normmd$Psd

  write.table(Amean[,],
              file=paste(fullOutputPath,"Amean.",FID,".txt",sep=""),
              sep="\t")
  write.table(Asd[,],
              file=paste(fullOutputPath,"Asd.",FID,".txt",sep=""),
              sep="\t")
  write.table(Pmean[,],
              file=paste(fullOutputPath,"Pmean.",FID,".txt",sep=""),
              sep="\t")
  write.table(Psd[,],
              file=paste(fullOutputPath,"Psd.",FID,".txt",sep=""),
              sep="\t")
    
    
  if (verbose) {
    message(paste("Computing chi2")) 
  }
    
  atmp <- as.matrix(Amean[,])
  ptmp <- as.matrix(Pmean[,])
  mtmp <- atmp %*% ptmp
  meanMock[,] = array(mtmp, c(nrow(D), ncol(D), 1))
  meanChi2 <- sum((D-meanMock[,])^2 / Sigma^2)

  # process the diagnostic file and add appropriate means
  ADiagFile <- paste("AResults",FID,".Diagnostics.txt",sep="")
  PDiagFile <- paste("PResults",FID,".Diagnostics.txt",sep="")

  ADiag <- read.table(ADiagFile, sep="\t", header=T, skip=10)
  PDiag <- read.table(PDiagFile, sep="\t", header=T, skip=10)

  chainKeepA  <- which(ADiag[,2]==0)
  meanOfChi2A <- mean(ADiag[chainKeepA,4])
  meanAtomsA  <- mean(ADiag[chainKeepA,5])

  chi2DiagA <- paste("<number of atoms>", meanAtomsA,
                     "<chi^2>:", meanOfChi2A,
                     "; chi^2 of mean:", meanChi2, sep=" ")
  cat("\n\nSummary\n", file=ADiagFile, append=TRUE)
  cat(chi2DiagA, file=ADiagFile, append=TRUE)
  cat("\n", file=ADiagFile, append=TRUE)
    
  chainKeepP  <- which(PDiag[,2]==0)
  meanOfChi2P <- mean(PDiag[chainKeepP,4])
  meanAtomsP  <- mean(PDiag[chainKeepP,5])

  chi2DiagP <- paste("<number of atoms>", meanAtomsP,
                     "<chi^2>:", meanOfChi2P,
                     "; chi^2 of mean:", meanChi2, sep=" ")
  cat("\n\nSummary\n", file=PDiagFile, append=TRUE)
  cat(chi2DiagP, file=PDiagFile, append=TRUE)
  cat("\n", file=PDiagFile, append=TRUE)
    
  if (fullOutputPath != "") {
     file.rename(ADiagFile, paste(fullOutputPath, basename(ADiagFile), sep=""))
     file.rename(PDiagFile, paste(fullOutputPath, basename(PDiagFile), sep=""))
  }

  # delete chain files or move to appropriate directory
  if (!keepChain) {
    message("Deleting chain files")
    file.remove(AFile)
    file.remove(PFile)
  } else {
    message(paste('Moving chain files to ', outputDir, sep=' '))
    if (fullOutputPath != "") {
       file.rename(AFile, paste(fullOutputPath, basename(AFile),sep=""))
       file.rename(PFile, paste(fullOutputPath, basename(PFile),sep=""))
    }
  }
  

  if (verbose) {
    message("Completed simulation.\n")
    message("Return list contains chi^2 of solution and normalized, mean and standard deviation of A and P matrices.\n")
  }

  # output list of data for each chain
  return(list(meanChi2=meanChi2,
              D=D, Sigma=Sigma,
              Amean=Amean, Asd=Asd,
              Pmean=Pmean, Psd=Psd,
              meanMock=meanMock))

}
                          
