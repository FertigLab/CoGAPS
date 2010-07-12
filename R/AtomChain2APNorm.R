# obtain mean and SD from output of data
AtomChain2APNorm <- function(Afile, Pfile,
                             nGene, nSample, nPattern,
                             verbose = TRUE) {

  if (verbose) {
    message(sprintf("Finding the distribution of A and P matrices from output files %s and %s\n.",
                  Afile, Pfile))
  }

  # determine the number of iterations from the file length

  nLinesA <- countLines(Afile)
  nLinesP <- countLines(Pfile)

  nIterA <- (nLinesA)
  nIterP <- (nLinesP)

  if ( (nIterA != nIterP) || nIterA <= 1 || nIterP <= 1 ) {
    stop(sprintf("Cannot normalize A and P matrices output at %d and %d iterations",
                 nIterA, nIterP))
  }

  if (verbose) {
    message("Computing mean of A and P matrices\n")
  }
  Anorm <- matrix(0, nrow = nGene,    ncol=nPattern)
  Pnorm <- matrix(0, nrow = nPattern, ncol=nSample)

  AConnection <- file(Afile,"r")
  if (!isOpen(AConnection, "r")) {
    stop(paste("Could not open",Afile,"for normalization.",sep=" "))
  }

  PConnection <- file(Pfile,"r")
  if (!isOpen(PConnection, "r")) {
    stop(paste("Could not open",Pfile,"for normalization.",sep=" "))
  }
  
  for (i in 1:nIterA) {
    AlineData <- scan(AConnection, sep="\t",nlines = 1, strip.white=T, quiet=T)
    PlineData <- scan(PConnection, sep="\t", nlines = 1, strip.white=T, quiet=T)

    if (length(AlineData) != nGene*nPattern) {
      stop(paste("Read invalid record for A matrix from", Afile, sep=" "))
    }

    if (length(PlineData) != nPattern*nSample) {
      stop(paste("Read invalid record for P matrix from", Pfile, sep=" "))
    }

    Atemp <- matrix(AlineData, nrow=nGene, ncol=nPattern)
    Ptemp <- matrix(PlineData, nrow=nPattern, ncol=nSample)

    normFactor <- apply(Ptemp,1,sum)
    Ptempnorm <- Ptemp / matrix(rep(normFactor, nSample),
                                nrow=nPattern, ncol=nSample)
    Atempnorm <- Atemp * matrix(rep(normFactor, nGene),
                                nrow=nGene, ncol=nPattern, byrow=T)

    
    Anorm <- Anorm + Atempnorm
    Pnorm <- Pnorm + Ptempnorm
  }
  Anorm <- Anorm / nIterA
  Pnorm <- Pnorm / nIterP

  close(AConnection)
  close(PConnection)
  
  if (verbose) {
    message("Computing sd of A and P matrices\n")
  }
  Asd <- matrix(0, nrow=nGene, ncol=nPattern)
  Psd <- matrix(0, nrow=nPattern, ncol=nSample)

  AConnection <- file(Afile,"r")
  if (!isOpen(AConnection, "r")) {
    stop(paste("Could not open",Afile,"for normalization.",sep=" "))
  }

  PConnection <- file(Pfile,"r")
  if (!isOpen(PConnection, "r")) {
    stop(paste("Could not open",Pfile,"for normalization.",sep=" "))
  }
  
  for (i in 1:nIterA) {
    AlineData <- scan(AConnection, sep="\t",nlines = 1, strip.white=T, quiet=T)
    PlineData <- scan(PConnection, sep="\t", nlines = 1, strip.white=T, quiet=T)

    if (length(AlineData) != nGene*nPattern) {
      stop(paste("Read invalid record for A matrix from", Afile, sep=" "))
    }

    if (length(PlineData) != nPattern*nSample) {
      stop(paste("Read invalid record for P matrix from", Pfile, sep=" "))
    }

    Atemp <- matrix(AlineData, nrow=nGene,    ncol=nPattern)
    Ptemp <- matrix(PlineData, nrow=nPattern, ncol=nSample)

    normFactor <- apply(Ptemp,1,sum)
    Ptempnorm <- Ptemp / matrix(rep(normFactor, nSample),
                                nrow=nPattern, ncol=nSample)
    Atempnorm <- Atemp * matrix(rep(normFactor, nGene),
                                nrow=nGene, ncol=nPattern, byrow=T)
    
    Asd <- Asd + (Atempnorm - Anorm)^2
    Psd <- Psd + (Ptempnorm - Pnorm)^2
  }
  Asd <- sqrt(Asd / (nIterA - 1))
  Psd <- sqrt(Psd / (nIterP - 1))
  
  close(AConnection)
  close(PConnection)

  return(list(Amean=Anorm, Pmean=Pnorm, Asd=Asd, Psd=Psd))
  
  
}
