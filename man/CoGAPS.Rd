\name{CoGAPS}
\alias{CoGAPS}
\title{CoGAPS driver script}

\description{
Runs the CoGAPS algorithm to infer underlying patterns in microarray data and their association to activity in gene sets.}

\usage{
  CoGAPS(data, unc, GStoGenes, outputDir, outputBase="", sep="\\t", 
         isPercentError=FALSE, numPatterns, MaxAtomsA=2^32, alphaA=0.01, 
         MaxAtomsP=2^32, alphaP=0.01, SAIter=1000000000, iter = 500000000, 
	 thin=-1, nPerm=500, verbose=TRUE, plot=FALSE, keepChain=FALSE)}

\arguments{
\item{data}{The matrix of m genes by n arrays of expression data.  The input can be either the data matrix itself or the file containing this data.  If the latter, CoGAPS will read in the data using \code{read.table(data, sep=sep, header=T, row.names=1)}.}
\item{unc}{The matrix of m genes by n arrays of uncertainty (standard deviation) for the expression data.  The input can be either a file containing the uncertainty (using the format from data), a matrix containing the uncertainty, or a constant value.  If unc is a constant value, it can represent either a constant uncertainty or a constant percentage of the values in data as determined by isPercentError.}
\item{GStoGenes}{List or data frame containing the genes in each gene set. If a list, gene set names are the list names and corresponding elements are the names of genes contained in each set. If a data frame, gene set names are in the first column and corresponding gene names are listed in rows beneath each gene set name.}
\item{numPatterns}{Number of patterns into which the data will be decomposed.  Must be less than the number of genes and number of arrays in the data.}
\item{outputDir}{Directory to which to output result and diagnostic files created by CoGAPS. (Use "" to output results to the current directory).}
\item{outputBase}{Prefix for all result and diagnostic files created by CoGAPS (optional; default="")}
\item{sep}{Text delimiter for tables in data and unc (if specified in file) and any output tables (optional; default="\\t")}
\item{isPercentError}{Boolean indicating whether constant value in unc is the value of the uncertainty or the percentage of the data that is the uncertainty.}
\item{MaxAtomsA}{Maximum number of atoms in the atomic domain used for the prior of the amplitude matrix in the decomposition (see Sibisi and Skilling, 1997).  The default value will typically be sufficient for most applications (optional; default=$2^32$).}
\item{alphaA}{Sparsity parameter reflecting the expected number of atoms per element of the amplitude matrix in the decomposition.  To enforce sparsity, this parameter should typically be less than one. (optional; default=0.01)}
\item{MaxAtomsP}{Maximum number of atoms in the atomic domain used for the prior of the pattern matrix in the decomposition (see Sibisi and Skilling, 1997).  The default value will typically be sufficient for most applications (optional; default=$2^32$).}
\item{alphaP}{Sparsity parameter reflecting the expected number of atoms per element of the pattern matrix in the decomposition.  To enforce sparsity, this parameter should typically be less than one. (optional; default=0.01)}
\item{SAIter}{Number of burn-in iterations for the MCMC matrix decomposition (optional; default=1000000000)}
\item{iter}{Number of iterations to represent the distribution of amplitude and pattern matrices with the MCMC matrix decomposition (optional; default=500000000)}
\item{thin}{Double whose integer part represents the number of iterations at which the samples are kept and decimal part provides an identifier for the output files from this implementation of CoGAPS.  If thin is an integer or not specified, this decimal file identifier is assigned randomly.  (optional; default=-1; code assigns number of iterations kept to be iter/1000 and file identifier to be runif(1)) }
\item{nPerm}{Number of permuations used for the null distribution in the gene set statistic. (optional; default=500)}
\item{verbose}{Boolean which specifies the amount of output to the user about the progress of the program. (optional; default=TRUE)}
\item{plot}{Boolean which specifies whether plots representing the resulting amplitude and pattern matrices should be made. (optional; default=FALSE)}
\item{keepChain}{Boolean which specifies if chain values of \eqn{{\bf{A}}} and \eqn{{\bf{P}}} are saved in outputDir (optional; default=FALSE).}
}

\details{
  CoGAPS first decomposes the data matrix using GAPS, \eqn{{\bf{D}}}, into a basis of underlying patterns and then determines the gene set activity in each of these patterns.

  The GAPS decomposition is achieved by finding amplitude and pattern matrices (\eqn{{\bf{A}}} and \eqn{{\bf{P}}}, respectively) for which \deqn{{\bf{D}} = {\bf{A}}{\bf{P}} + \Sigma,} where \eqn{\Sigma} is the matrix of uncertainties given by unc.  The matrices \eqn{\bf{A}} and \eqn{\bf{P}} are assumed to have the atomic prior described in Sibisi and Skilling (1997) and are found with MCMC sampling implemented within JAGS.

  Then, the patterns identified in the columns of \eqn{\bf{P}} are linked to activity in each of the gene sets specified in GStoGenes using a novel z-score based statistic developed in Ochs et al. (2009).  Specifically, the z-score for pattern \eqn{p} and gene set \eqn{G_{i}} containing $G$ total genes is given by \deqn{Z_{i,p} = \frac{1}{G} \sum_{g in \mathcal{G_{i}}} {\frac{{\bf{A}_{gp}}}{Asd_{gp}}},}  
where \eqn{g} indexes the genes in the set and \eqn{Asd_{gp}} is the standard deviation of \eqn{{\bf{A}}_{gp}} obtained from MCMC sampling.  CoGAPS then uses the specified \code{nPerm} random sample tests to compute a consistent p value estimate from that z score.  Note that the data from Ochs et al. (2009) are provided with this package in GIST_TS_20084.RData and TFGSList.RData are also provided with this package for further validation with nIter=5e+07.
}

\value{
  A list containing:
  \item{D}{Microarray data matrix.}
  \item{Sigma}{Data matrix with uncertainty of D.}
  \item{Amean}{Sampled mean value of the amplitude matrix \eqn{{\bf{A}}}.}
  \item{Asd}{Sampled standard deviation of the amplitude matrix \eqn{{\bf{A}}}.}
  \item{Pmean}{Sampled mean value of the pattern matrix \eqn{\bf{P}}.}
  \item{Psd}{Sampled standard deviation of the pattern matrix \eqn{\bf{P}}.}
  \item{meanMock}{Mock data obtained from matrix decomposition for sampled mean values (= Amean \%*\% Pmean).}
  \item{meanChi2}{\eqn{\chi^2} value for the sampled mean values (Amean and Pmean) of the matrix decomposition.}
  \item{GSUpreg}{p-values for upregulation of each gene set in each pattern.}
  \item{GSDownreg}{p-values for downregulation of each gene set in each pattern.}
  \item{GSActEst}{p-values for activity of each gene set in each pattern.}
}

\note{Running GAPS will create the folder ouptutDir, create diagnostic files with \eqn{\chi^{2}} and number of atoms, files with the mean and standard deviation of \eqn{{\bf{A}}} and \eqn{{\bf{P}}}, files with p-values for upregulation/downregulation/activity of each gene set, and optionally values of \eqn{{\bf{A}}} and \eqn{{\bf{P}}} from the MCMC chain.}

\author{Elana J. Fertig \email{ejfertig@jhmi.edu}}

\references{
M.F. Ochs, L. Rink, C. Tarn, S. Mburu, T. Taguchi, B. Eisenberg, and A.K. Godwin.  (2009) Detection and treatment-induced changes in signaling pathways in gastrointestinal stromal tumors using transcriptomic data.  Cancer Research, 69:9125-9132.

M. Plummer. (2003) JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling. In K. Hornik, F. Leisch, and A. Zeileis, editors, Proceedings of the Third International Workshop on Distributed Statistical Computing, Vienna, Austria.

S. Sibisi and J. Skilling. (1997) Prior distributions on measure space.  Journal of the Royal Statistical Society, B, 59:217-235.
}

\examples{
\dontrun{
## Load data
data(EasySimGS)

## Run the CoGAPS matrix decomposition
nIter <- 5e+05 
results <- CoGAPS(data=DGS, unc=0.01, isPercentError=FALSE,
                  GStoGenes=gs,
                  numPatterns=3,
                  SAIter = 2*nIter, iter = nIter,
                  outputDir='GSResults', plot=TRUE)
}
}

\seealso{\code{\link{GAPS}}, \code{\link{calcCoGAPSStat}}}
\keyword{misc}
