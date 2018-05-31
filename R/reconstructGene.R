#' Reconstruct Gene
#'
#' @param A A matrix estimates
#' @param P corresponding P matrix estimates
#' @param genes an index of the gene or genes of interest
#' @return the D' estimate of a gene or set of genes
#' @examples
#' data(SimpSim)
#' reconstructGene(SimpSim.result$Amean, SimpSim.result$Pmean)
#' @export
reconstructGene<-function(A, P, genes=NA)
{
    Dneu <- A %*% P
    if (!is.na(genes))
        Dneu <- Dneu[genes,]
    return(Dneu)
}
