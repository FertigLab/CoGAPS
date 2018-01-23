# Output: list with gene set statistics
#' Calculate Gene Set Statistics
#'
#' @details calculates the gene set statistics for each
#'  column of A using a Z-score from the elements of the A matrix,
#'  the input gene set, and permutation tests
#' @param Amean A matrix mean values
#' @param Asd A matrix standard deviations
#' @param GStoGenes data.frame or list with gene sets
#' @param numPerm number of permutations for null
#' @export
calcCoGAPSStat <- function (Amean, Asd, GStoGenes, numPerm=500)
{
    # test for std dev of zero, possible mostly in simple simulations
    if (sum(Asd==0) > 0)
        Asd[Asd==0] <- .Machine$double.eps

    # calculate Z scores
    zMatrix <- calcZ(Amean,Asd)

    # check input arguments
    if (!is(GStoGenes, "data.frame") && !is(GStoGenes, "list") && !is(GStoGenes,"GSA.genesets"))
        stop("GStoGenes must be a data.frame,GSA.genesets, or list with format specified in the users manual.")

    if (is(GStoGenes, "GSA.genesets"))
    {
        names(GStoGenes$genesets) <- GStoGenes$geneset.names
        GStoGenes <- GStoGenes$genesets
    }

    if (is(GStoGenes, "list"))
        GStoGenesList <- GStoGenes
    else
        GStoGenesList <- list()

    for (i in 1:dim(GStoGenes)[2])
    {
        GStoGenesList[[as.character(colnames(GStoGenes)[i])]] <- as.character(unique(GStoGenes[,i]))
    }

    # get dimensions
    numGS   <- length(names(GStoGenesList))
    numPatt <- dim(zMatrix)[2]
    numG    <- dim(zMatrix)[1]+0.9999

    # initialize matrices
    statsUp   <- matrix(ncol = numGS, nrow = numPatt)
    statsDown <- matrix(ncol = numGS, nrow = numPatt)
    actEst    <- matrix(ncol = numGS, nrow = numPatt)
    results   <- list()
    zPerm     <- matrix(ncol=numPerm,nrow=numPatt)

    # do permutation test
    for (gs in 1:numGS)
    {
        genes <- GStoGenesList[[names(GStoGenesList)[gs]]]
        index <- rownames(zMatrix) %in% genes
        zValues <- zMatrix[index,1]
        numGenes <- length(zValues)
        label <- as.character(numGenes)

        if (!any(names(results)==label))
        {
            for (p in 1:numPatt)
            {
                for (j in 1:numPerm)
                {
                    temp <- floor(runif(numGenes,1,numG))
                    temp2 <- zMatrix[temp,p]
                    zPerm[p,j] <- mean(temp2)
                }
            }
            results[[label]] <- zPerm
        }
    }

    # get p-value
    for (p in 1:numPatt)
    {
        for (gs in 1:numGS)
        {
            genes <- GStoGenesList[[names(GStoGenesList)[gs]]]
            index <- rownames(zMatrix) %in% genes
            zValues <- zMatrix[index,p]
            zScore <- mean(zValues)

            numGenes <- length(zValues)
            label <- as.character(numGenes)

            permzValues <- results[[label]][p,]
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

    results[['GSUpreg']] <- statsUp
    results[['GSDownreg']] <- statsDown
    results[['GSActEst']] <- actEst

    return(results)
}

# calcGeneGSStat: calculate probability gene belongs in gene set
# History: v1.0 EJF original CoGAPS

# Inputs: Amean - A matrix mean values
#         Asd - A matrix standard deviations
#         GStoGenes - data.frame or list with gene sets
#         numPerm - number of permutations for null
#         Pw - weighting on genes
#         nullGenes - adjust genes

# Output: environment with gene set statistics
#       NOTE: should make into list object, env historical

#'\code{calcGeneGSStat} calculates the probability that a gene
#'listed in a gene set behaves like other genes in the set within
#'the given data set
#'
#'@param Amean A matrix mean values
#'@param Asd A matrix standard deviations
#'@param GSGenes data.frame or list with gene sets
#'@param numPerm number of permutations for null
#'@param Pw weight on genes
#'@param nullGenes - logical indicating gene adjustment
#'@export
calcGeneGSStat  <- function(Amean, Asd, GSGenes, numPerm, Pw=rep(1,ncol(Amean)), nullGenes=FALSE)
{
    gsStat <- calcCoGAPSStat(Amean, Asd, data.frame(GSGenes), numPerm=numPerm)
    gsStat <- gsStat$GSUpreg

    gsStat <- -log(gsStat)

    if (!all(is.na(Pw)))
    {
        if (length(Pw) != length(gsStat))
        {
            stop('Invalid weighting')
        }
        gsStat <- gsStat*Pw
    }

    if (nullGenes)
    {
        ZD <- Amean[setdiff(row.names(Amean), GSGenes),] /
            Asd[setdiff(row.names(Amean), GSGenes),]
    }
    else
    {
        ZD <- Amean[GSGenes,]/Asd[GSGenes,]
    }
    outStats <- apply(sweep(ZD,2,gsStat,FUN="*"),1,sum) / (sum(gsStat))
    outStats <- outStats / apply(ZD,1,sum)
    outStats[which(apply(ZD,1,sum) < .Machine$double.eps)] <- 0

    if (sum(gsStat) < .Machine$double.eps)
        return(0)
    return(outStats)
}

computeGeneGSProb <- function(Amean, Asd, GSGenes, Pw=rep(1,ncol(Amean)),
numPerm=500, PwNull=FALSE)
{
    geneGSStat <- calcGeneGSStat(Amean=Amean, Asd=Asd, Pw=Pw,
        GSGenes=GSGenes, numPerm=numPerm)

    if (PwNull)
    {
        permGSStat <- calcGeneGSStat(Amean=Amean, Asd=Asd, GSGenes=GSGenes,
            numPerm=numPerm, Pw=Pw, nullGenes=T)
    }
    else
    {
        permGSStat <- calcGeneGSStat(Amean=Amean, Asd=Asd, GSGenes=GSGenes,
        numPerm=numPerm, nullGenes=T)
    }

    return(sapply(GSGenes, function(x)
        length(which(permGSStat > geneGSStat[x])) / length(permGSStat)))
}

#' Compute Z-Score Matrix
#'
#' @details calculates the Z-score for each element based on input mean
#'  and standard deviation matrices
#' @param meanMat matrix of mean values
#' @param sdMat matrix of standard deviation values
#' @export
calcZ <- function (meanMat, sdMat)
{
    if (nrow(meanMat) != nrow(sdMat))
        stop("Number of rows in the mean and standard deviation of A do not agree.")
    if (ncol(meanMat) != ncol(sdMat))
        stop("Number of columns in the mean and standard deviation of A do not agree.")

    zMat <- meanMat / sdMat
    rownames(zMat) <- rownames(meanMat)
    colnames(zMat) <- colnames(meanMat)

    return(zMat)
}

#' Reconstruct Gene
#'
#' @param A A matrix estimates
#' @param P corresponding P matrix estimates
#' @param genes an index of the gene or genes of interest
#' @return the D' estimate of a gene or set of genes
#' @export
reconstructGene<-function(A=NA, P=NA, genes=NA)
{
    Dneu <- A %*% P
    if (!is.na(genes))
        Dneu <- Dneu[genes,]
    return(Dneu)
}

