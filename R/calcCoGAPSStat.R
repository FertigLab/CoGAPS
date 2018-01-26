#' Calculate Gene Set Statistics
#'
#' @details calculates the gene set statistics for each
#' column of A using a Z-score from the elements of the A matrix,
#' the input gene set, and permutation tests
#' @param Amean A matrix mean values
#' @param Asd A matrix standard deviations
#' @param GStoGenes data.frame or list with gene sets
#' @param numPerm number of permutations for null
#' @return gene set statistics for each column of A
#' @examples
#' # Load the sample data from CoGAPS
#' data(SimpSim)
#' # Run calcCoGAPSStat with the correct arguments from 'results'
#' calcCoGAPSStat(SimpSim.result$Amean, SimpSim.result$Asd,
#' GStoGenes=GSets, numPerm=500)
#' @export
calcCoGAPSStat <- function (Amean, Asd, GStoGenes, numPerm=500)
{
    # test for std dev of zero, possible mostly in simple simulations
    if (sum(Asd==0) > 0)
        Asd[Asd==0] <- 1e-6

    # calculate Z scores
    zMatrix <- calcZ(Amean,Asd)

    if (is(GStoGenes, "GSA.genesets"))
    {
        names(GStoGenes$genesets) <- GStoGenes$geneset.names
        GStoGenes <- GStoGenes$genesets
    }
    else if (is(GStoGenes, "list"))
    {
        GStoGenesList <- GStoGenes
    }
    else if (is(GStoGenes, "data.frame"))
    {
        GStoGenesList <- list()
        for (i in 1:dim(GStoGenes)[2])
        {
            GStoGenesList[[as.character(colnames(GStoGenes)[i])]] <- 
                as.character(unique(GStoGenes[,i]))
        }
    }
    else
    {
        stop(paste("GStoGenes must be a data.frame, GSA.genesets, or list with",
            "format specified in the users manual."))
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

