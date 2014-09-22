\name{computeGeneGSProb}
\alias{computeGeneGSProb}
\alias{geneGSProb}
\title{CoGAPS gene membership statistic}

\description{Computes the p-value for gene set membership using the CoGAPS-based statistics developed in Fertig et al. (2012).  This statistic refines set membership for each candidate gene in a set specified in \code{GSGenes} by comparing the inferred activity of that gene to the average activity of the set.  Specifically, we compute the following summary statistic for each gene \eqn{g} that is a candidate member of gene set \eqn{G}:
\deqn{S_{g,G} = (\sum_{p} -log(Pr_{G,p})Pw[p](A_{gp}/\sigma_{gp})) / \sum_{p}-log(Pr_{G,p})Pw[p],}
where \eqn{p} indexes each of the patterns, \eqn{Pr_{G,p}} is the probability that gene set \eqn{G} is upregulated computed with \code{\link{calcCoGAPSStat}}, \eqn{A_{gp}} is the mean amplitude matrix from the \code{GAPS} matrix factorization, Pw[p] is a prior weighting for each pattern based upon the context to which that pattern relates, and \eqn{\sigma_{gp}} is the standard deviation of the amplitude matrix. P-values are formulated from a permutation test comparing the value of \eqn{S_{g,G}} for genes in \code{GSGenes} relative to the value of \eqn{S_{g,G}} \code{numPerm} random gene sets with the same number of targets.
}

\usage{computeGeneGSProb(Amean, Asd, GStoGenes, Pw=rep(1,ncol(Amean)),numPerm=500,PwNull=F)}

\arguments{
\item{Amean}{Sampled mean value of the amplitude matrix \eqn{{\bf{A}}}.  \code{row.names(Amean)} must correspond to the gene names contained in GStoGenes.}
\item{Asd}{Sampled standard deviation of the amplitude matrix \eqn{{\bf{A}}}.}
\item{GStoGenes}{Vector containing the prior estimate of members of the gene set of interest.}
\item{Pw}{Vector containing the weight to assign each pattern in the gene statistic assumed to be computed from the association of the pattern with samples in a given context (optional: default=1 giving all patterns equal weight).}
\item{numPerm}{Number of permuations used for the null distribution in the gene set statistic. (optional; default=500)}
\item{PwNull}{Logical value. If TRUE, use pattern weighting in Pw when computing the null distribution for the statistic.  If FALSE, do not use the pattern weighting so that the null is context independent. (optional; default=F)}
}

\value{
  A vector of length GStoGenes containing the p-values of set membership for each gene containined in the set specified in GStoGenes.
}

\examples{
\dontrun{

#################################################
# Results for GIST data in Fertig et al. (2012) #
#################################################

# load the data
data('GIST_TS_20084')
data('TFGSList')

# define transcription factors of interest based on Ochs et al. (2009)
TFs <- c("c.Jun", 'NF.kappaB', 'Smad4', "STAT3", "Elk.1", "c.Myc", "E2F.1", 
         "AP.1", "CREB", "FOXO", "p53", "Sp1")

# run the GAPS matrix factorization
nIter <- 10000
results <- CoGAPS(GIST.D, GIST.S, tf2ugFC,
                  nFactor=5,
                  nEquil=nIter, nSample=nIter,
                  plot=FALSE)

# set membership statistics
permTFStats <- list()
for (tf in TFs) {
     genes <- levels(tf2ugFC[,tf])
     genes <- genes[2:length(genes)]
     permTFStats[[tf]] <- computeGeneTFProb(Amean = GISTResults$Amean, 
                                            Asd = GistResults$Asd, genes) 
}

}
}

\author{Elana J. Fertig \email{ejfertig@jhmi.edu}}

\references{
E.J. Fertig, A.V. Favorov, and M.F. Ochs (2013) Identifying context-specific transcription factor targets from prior knowledge and gene expression data. 2012 IEEE Nanobiosciences. 
}

\seealso{\code{\link{calcCoGAPSStat}}}
\keyword{misc}
