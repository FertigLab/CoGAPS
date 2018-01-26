#' Probability Gene Belongs in Gene Set
#'
#' @details calculates the probability that a gene
#' listed in a gene set behaves like other genes in the set within
#' the given data set
#' @param Amean A matrix mean values
#' @param Asd A matrix standard deviations
#' @param GSGenes data.frame or list with gene sets
#' @param numPerm number of permutations for null
#' @param Pw weight on genes
#' @param nullGenes logical indicating gene adjustment
#' @return gene similiarity statistic
#' @examples
#' # Load the sample data from CoGAPS
#' data('SimpSim')
#' # Run calcGeneGSStat with the correct arguments from 'results'
#' calcGeneGSStat(SimpSim.result$Amean, SimpSim.result$Asd, 
#' GSGenes=GSets[[1]], numPerm=500)
#' @export
calcGeneGSStat  <- function(Amean, Asd, GSGenes, numPerm, Pw=rep(1,ncol(Amean)),
nullGenes=FALSE)
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
    {
        return(0)
    }
    return(outStats)
}

#' Compute Gene Probability
#'
#' @details Computes the p-value for gene set membership using the CoGAPS-based
#' statistics developed in Fertig et al. (2012).  This statistic refines set
#' membership for each candidate gene in a set specified in \code{GSGenes} by
#' comparing the inferred activity of that gene to the average activity of the
#' set.
#' @param Amean A matrix mean values
#' @param Asd A matrix standard deviations
#' @param GSGenes data.frame or list with gene sets
#' @param Pw weight on genes
#' @param numPerm number of permutations for null
#' @param PwNull - logical indicating gene adjustment
#' @return A vector of length GSGenes containing the p-values of set membership
#' for each gene containined in the set specified in GSGenes.
#' @examples
#' # Load the sample data from CoGAPS
#' data('SimpSim')
#' # Run calcGeneGSStat with the correct arguments from 'results'
#' calcGeneGSStat(SimpSim.result$Amean, SimpSim.result$Asd, 
#' GSGenes=GSets[[1]], numPerm=500)
#' @export
computeGeneGSProb <- function(Amean, Asd, GSGenes, Pw=rep(1,ncol(Amean)),
numPerm=500, PwNull=FALSE)
{
    geneGSStat <- calcGeneGSStat(Amean=Amean, Asd=Asd, Pw=Pw,
        GSGenes=GSGenes, numPerm=numPerm)

    if (PwNull)
    {
        permGSStat <- calcGeneGSStat(Amean=Amean, Asd=Asd,
            GSGenes=GSGenes, numPerm=numPerm, Pw=Pw,
            nullGenes=TRUE)
    }
    else
    {
        permGSStat <- calcGeneGSStat(Amean=Amean, Asd=Asd,
            GSGenes=GSGenes, numPerm=numPerm,
            nullGenes=TRUE)
    }

    finalStats <- sapply(GSGenes, function(x)
        length(which(permGSStat > geneGSStat[x])) / length(permGSStat))

    return(finalStats)
}
