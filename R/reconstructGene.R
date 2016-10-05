#' reconstruct Gene
#'
#' @param A A matrix estimates
#' @param P corresponding P matrix estimates
#' @param genes an indx of the gene or genes of interest. If \code{NA}, the default, all genes contained in A will be returned.
#'
#' @return the D' estimate of a gene or set of genes
#' @export
#'
reconstructGene<-function(A=AP$Amean,P=AP$Pmean,genes=NA){
  Dneu <- A%*%P
  if(!is.na(genes)){Dneu <- Dneu[genes,]}
  return(Dneu)
}

