\name{computeGeneGSProb}
\alias{computeGeneGSProb}
\alias{geneGSProb}
\title{CoGAPS gene membership statistic}

\description{Computes the p-value for gene set membership based on the CoGAPS inferred activity of genes relative to the activity of prior estimates for members of the gene sets.}

\usage{computeGeneGSProb(Amean, Asd, GSGenes, numPerm=500)}

\arguments{
\item{Amean}{Sampled mean value of the amplitude matrix \eqn{{\bf{A}}}.  \code{row.names(Amean)} must correspond to the gene names contained in GStoGenes.}
\item{Asd}{Sampled standard deviation of the amplitude matrix \eqn{{\bf{A}}}.}
\item{GSGenes}{Vector containing the prior estimate of members of the gene set of interest.}
\item{numPerm}{Number of permuations used for the null distribution in the gene set statistic. (optional; default=500)}
}

\value{
  A vector of length GSGenes containing the p-values of gene set membership.
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
