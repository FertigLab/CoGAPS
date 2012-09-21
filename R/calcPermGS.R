calcGeneGSStat  <- function(Amean, Asd, GSGenes, numPerm, nullGenes=F) {

   # Compute the gene set statistics of upregulation
   gsStat <- calcCoGAPSStat(Amean, Asd, data.frame(GSGenes), 
  			    numPerm=numPerm)
   gsStats <- -log(gsStat$GSUpreg)

   if (sum(gsStat) < .Machine$double.eps) {
      return(0.)
   }

   # get the ranked statistic for targets of the gene set		
   if (nullGenes) {
      # genes outside of the set serve as a null distribution
      ZD <- Amean[setdiff(row.names(Amean), GSGenes),] /
            Asd[setdiff(row.names(Amean), GSGenes),]
   } else {
      # statistics for genes within the set
      ZD <- Amean[GSGenes,]/Asd[GSGenes,]
   }
  
   # compute the statistic
   outStats <- apply(sweep(ZD,2,gsStat,FUN="*"),1,sum) / (sum(gsStat))
   outStats <- outStats / apply(ZD,1,sum) # normalization
   outStats[which(apply(ZD,1,sum) < .Machine$double.eps)] <- 0.
	
   return(outStats)
	
}


computeGeneGSProb <- function(Amean, Asd, GSGenes, numPerm=500) {

   # get the gene membership statistic for genes within the set
   geneGSStat <- calcGeneGSStat(Amean=Amean, Asd=Asd,
				GSGenes=GSGenes, numPerm=numPerm)
	
   # get the gene membership statistic for null genes
   nullGSStat <- calcGeneGSStat(Amean=Amean, Asd=Asd,
			        GSGenes=GSGenes, numPerm=numPerm,
				nullGenes=T)
	
   # perform the permutation test to get the p-values
   finalStats <- sapply(GSGenes,
			function(x){length(which(nullGSStat > geneGSStat[x])) /
                                    length(nullGSStat)})
	
   return(finalStats)
}
