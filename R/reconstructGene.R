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

