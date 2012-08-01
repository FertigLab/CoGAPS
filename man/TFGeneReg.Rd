\name{TFGeneReg}
\docType{data}
\alias{TFGeneReg}
\alias{TFSimData}
\alias{GSSimData}
\title{Simulated dataset to quantify gene set membership.}
\description{Simulated data and components used to generate it resulting from the differential activity of four simulated gene sets (\code{TFGeneReg$TFGeneReg}) in different samples (\code{TFGeneReg$P}).}
\usage{TFGeneReg}
\format{A \code{\link{list}} containing:
	  \item{A}{Matrix of 100 rows and 4 columns representing the simulated amplitude matrix for activity of each of the four simulated gene sets in each pattern.}
	  \item{D}{Matrix of 100 rows and 20 columns containing simulated data generated with \code{M + 0.1*pmax(TFGeneReg$M,1)*matrix(rnorm(length(TFGeneReg$M)),nrow=nrow(TFGeneReg$M))}.}
	  \item{M}{Matrix of 100 rows and 20 columns containing noise-free simulated data generated with \code{TFGeneReg$A %*% TFGeneReg$P}.}
	  \item{P}{Matrix of 4 columns and 20 columns representing relative activity of each of the four gene sets in \code{TFGeneReg$TFGeneReg} in each of the 20 samples.}
	  \item{TFGeneReg}{List containing genes and relative activity for each of four gene sets used to formulate the amplitude matrix \code{TFGeneReg$A}.}
}
\references{EJ Fertig, AV Favorov, and Ochs MF (2012) Identifying context-specific transcription factor targets from prior knowledge and gene expression data. 2012 IEEE International Conference on Bioinformatics and Biomedicine.}
\keyword{datasets}
