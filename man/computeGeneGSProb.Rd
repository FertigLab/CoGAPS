\name{computeGeneGSProb}
\alias{computeGeneGSProb}
\alias{geneGSProb}
\title{CoGAPS gene membership statistic}

\description{Computes the p-value for gene set membership using the CoGAPS-based statistics developed in Fertig et al. (2012).  This statistic refines set membership for each candidate gene in a set specified in \code{GSGenes} by comparing the inferred activity of that gene to the average activity of the set.  Specifically, we compute the following summary statistic for each gene \eqn{g} that is a candidate member of gene set \eqn{G}:
\deqn{S_{g,G} = (\sum_{p} -log(Pr_{G,p})(A_{gp}/\sigma_{gp})) / \sum_{p}-log(Pr_{G,p}),}
where \eqn{p} indexes each of the patterns, \eqn{Pr_{G,p}} is the probability that gene set \eqn{G} is upregulated computed with \code{\link{calcCoGAPSStat}}, \eqn{A_{gp}} is the mean amplitude matrix from the \code{GAPS} matrix factorization, and \eqn{\sigma_{gp}} is the standard deviation of the amplitude matrix. P-values are formulated from a permutation test comparing the value of \eqn{S_{g,G}} for genes in \code{GSGenes} relative to the value of \eqn{S_{g,G}} \code{numPerm} random gene sets with the same number of targets.
}

\usage{computeGeneGSProb(Amean, Asd, GSGenes, numPerm=500)}

\arguments{
\item{Amean}{Sampled mean value of the amplitude matrix \eqn{{\bf{A}}}.  \code{row.names(Amean)} must correspond to the gene names contained in GStoGenes.}
\item{Asd}{Sampled standard deviation of the amplitude matrix \eqn{{\bf{A}}}.}
\item{GSGenes}{Vector containing the prior estimate of members of the gene set of interest.}
\item{numPerm}{Number of permuations used for the null distribution in the gene set statistic. (optional; default=500)}
}

\value{
  A vector of length GSGenes containing the p-values of set membership for each gene containined in the set specified in GSGenes.
}

\examples{
\dontrun{
##########################################
# Simulated data in Fertig et al. (2012) #
##########################################

## Load data
data(TFSimData)

## Run the CoGAPS matrix decomposition
nIter <- 5e+07 

results <- GAPS(data=TFGeneReg$D, 
                unc=0.1*pmax(TFGeneReg$M,1), 
                isPercentError=FALSE,
                numPatterns=4,
                SAIter = 2*nIter, iter = nIter,
                outputDir='GSResults')

# compute the probablity of membership of each gene in each set
TFtargets <- lapply(TFGeneReg$TFGeneReg,names)
TFGenesP <- lapply(TFtargets, function(x){
   computeGeneGSProb(Amean=results$Amean, Asd=results$Asd, GSGenes=x)
})

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
nIter <- 5e7

GISTResults <- GAPS(data=GIST.D, unc=GIST.S, 
                    numPatterns=5, outputDir = 'GISTGSCoGAPS',
                    isPercentError=F, SAIter=2*nIter, iter=nIter)

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
E.J. Fertig, A.V. Favorov, and M.F. Ochs (2012) Identifying context-specific transcription factor targets from prior knowledge and gene expression data. 2012 IEEE Conference on Bioinformatics and Biomedicine. 
}

\seealso{\code{\link{calcCoGAPSStat}}}
\keyword{misc}
