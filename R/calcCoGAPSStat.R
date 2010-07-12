calcCoGAPSStat <- function (Amean, Asd, GStoGenes, numPerm=500) {

  # calculate Z scores
  zMatrix <- calcZ(Amean,Asd)
  
  # compute the p-value for each gene set belonging to each pattern

  # check input arguments
  if (!is(GStoGenes, "data.frame") && !is(GStoGenes, "list")) {
    stop("GStoGenes must be a data.frame or list with format specified in the users manual.")
  }
	
  if (is(GStoGenes, "list")) {
    GStoGenesList <- GStoGenes
  } else {
    GStoGenesList <- list()
    for (i in 1:dim(GStoGenes)[2]) {
      GStoGenesList[[as.character(colnames(GStoGenes)[i])]] <- as.character(unique(GStoGenes[,i]))
    }
  }
  
  # get dimensions
  numGS   <- length(names(GStoGenesList))
  numPatt <- dim(zMatrix)[2]
  numG    <- dim(zMatrix)[1]+0.9999

  # initialize matrices
  statsUp   <- matrix(ncol = numGS, nrow = numPatt)
  statsDown <- matrix(ncol = numGS, nrow = numPatt)
  actEst    <- matrix(ncol = numGS, nrow = numPatt)
  results   <- new.env()
  zPerm     <- matrix(ncol=numPerm,nrow=numPatt)

  # do permutation test
  for (gs in 1:numGS) {
    genes <- GStoGenesList[[names(GStoGenesList)[gs]]]
    index <- rownames(zMatrix) %in% genes
    zValues <- zMatrix[index,1]
    numGenes <- length(zValues)
    label <- as.character(numGenes)
    
    if (!exists(label,env=results)) {
      for (p in 1:numPatt) {
        for (j in 1:numPerm) {
          temp <- floor(runif(numGenes,1,numG))
          temp2 <- zMatrix[temp,p]
          zPerm[p,j] <- mean(temp2)
        }
      }
      assign(label,zPerm,env=results)
    }
  }

  # get p-value
  for (p in 1:numPatt) {

    for (gs in 1:numGS) {

      genes <- GStoGenesList[[names(GStoGenesList)[gs]]]
      index <- rownames(zMatrix) %in% genes
      zValues <- zMatrix[index,p]
      zScore <- mean(zValues)

      numGenes <- length(zValues)
      label <- as.character(numGenes)

      permzValues <- get(label, env=results)[p,]
      ordering <- order(permzValues)
      permzValues <- permzValues[ordering]

      statistic <- sum(zScore > permzValues)
      statUpReg <- 1 - statistic/length(permzValues) 
      statsUp[p,gs] <- max(statUpReg, 1/numPerm)

      statistic <- sum(zScore < permzValues)
      statDownReg <- 1 - statistic/length(permzValues) 
      statsDown[p,gs] <- max(statDownReg,1/numPerm)

      activity <- 1 - 2*max(statUpReg, 1/numPerm)
      actEst[p,gs] <- activity
    }

  }

  # format output
  colnames(statsUp) <- names(GStoGenesList)
  colnames(statsDown) <- names(GStoGenesList)
  colnames(actEst) <- names(GStoGenesList)

  rownames(statsUp) <- colnames(zMatrix)
  rownames(statsDown) <- colnames(zMatrix)
  rownames(actEst) <- colnames(zMatrix)

  assign('GSUpreg', statsUp, env=results)
  assign('GSDownreg', statsDown, env=results)
  assign('GSActEst', actEst, env=results)
  
  return(results)
}
