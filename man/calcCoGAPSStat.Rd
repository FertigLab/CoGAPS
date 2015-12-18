\name{calcCoGAPSStat}
\alias{calcCoGAPSStat}
\title{CoGAPS gene set statistic}

\description{
Computes the p-value for the association of underlying patterns from microarray data to activity in gene sets.}

\usage{
  calcCoGAPSStat(Amean, Asd, GStoGenes, numPerm=500)}

\arguments{
\item{Amean}{Sampled mean value of the amplitude matrix \eqn{{{A}}}.  \code{row.names(Amean)} must correspond to the gene names contained in GStoGenes.}
\item{Asd}{Sampled standard deviation of the amplitude matrix \eqn{{{A}}}.}
\item{GStoGenes}{List or data frame containing the genes in each gene set. If a list, gene set names are the list names and corresponding elements are the names of genes contained in each set. If a data frame, gene set names are in the first column and corresponding gene names are listed in rows beneath each gene set name.}
\item{numPerm}{Number of permuations used for the null distribution in the gene set statistic. (optional; default=500)}
}

\details{
  This script links the patterns identified in the columns of \eqn{{P}}  to activity in each of the gene sets specified in GStoGenes using a novel z-score based statistic developed in Ochs et al. (2009).  Specifically, the z-score for pattern \eqn{p} and gene set \eqn{G_{i}} containing \eqn{G} total genes is given by \deqn{Z_{i,p} = \frac{1}{G} \sum_{g in G_{i}}A_{gp} / \sigma_{gp}}, where \eqn{g} indexes the genes in the set and \eqn{\sigma_{gp}} is the standard deviation of \eqn{{{A}}_{gp}} obtained from MCMC sampling.  CoGAPS then uses the specified \code{numPerm} random sample tests to compute a consistent p value estimate from that z score.
}

\value{
  A list containing:
  \item{GSUpreg}{p-values for upregulation of each gene set in each pattern.}
  \item{GSDownreg}{p-values for downregulation of each gene set in each pattern.}
  \item{GSActEst}{p-values for activity of each gene set in each pattern.}
}

\author{Elana J. Fertig \email{ejfertig@jhmi.edu}}

\references{
M.F. Ochs, L. Rink, C. Tarn, S. Mburu, T. Taguchi, B. Eisenberg, and A.K. Godwin.  (2009) Detection and treatment-induced changes in signaling pathways in gastrointestinal stromal tumors using transcriptomic data.  Cancer Research, 69:9125-9132.
}

\seealso{\code{\link{CoGAPS}}}
\keyword{misc}
