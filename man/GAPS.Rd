\name{GAPS}
\alias{GAPS}
\title{GAPS matrix decomposition script}

\description{Decomposes microarray data into underlying patterns and corresponding amplitude.}

\usage{
  GAPS(data, unc, outputDir, outputBase="", sep="\\t", isPercentError=FALSE, 
       numPatterns, MaxAtomsA=2^32, alphaA=0.01, MaxAtomsP=2^32, alphaP=0.01, 
       SAIter=1000000000, iter = 500000000, thin=-1, verbose=TRUE,
       keepChain=FALSE)}

\arguments{
\item{data}{The matrix of m genes by n arrays of expression data.  The input can be either the data matrix itself or the file containing this data.  If the latter, GAPS will read in the data using \code{read.table(data, sep=sep, header=T, row.names=1)}.}
\item{unc}{The matrix of m genes by n arrays of uncertainty (standard deviation) for the expression data.  The input can be either a file containing the uncertainty (using the format from data), a matrix containing the uncertainty, or a constant value.  If unc is a constant value, it can represent either a constant uncertainty or a constant percentage of the values in data as determined by isPercentError.}
\item{numPatterns}{Number of patterns into which the data will be decomposed.  Must be less than the number of genes and number of arrays in the data.}
\item{outputDir}{Directory to which to output result and diagnostic files created by GAPS. (Use "" to output results to the current directory).}
\item{outputBase}{Prefix for all result and diagnostic files created by GAPS (optional; default="")}
\item{sep}{Text delimiter for tables in data and unc (if specified in file) and any output tables (optional; default="\\t")}
\item{isPercentError}{Boolean indicating whether constant value in unc is the value of the uncertainty or the percentage of the data that is the uncertainty.}
\item{MaxAtomsA}{Maximum number of atoms in the atomic domain used for the prior of the amplitude matrix in the decomposition (see Sibisi and Skilling, 1997).  The default value will typically be sufficient for most applications (optional; default=$2^32$).}
\item{alphaA}{Sparsity parameter reflecting the expected number of atoms per element of the amplitude matrix in the decomposition.  To enforce sparsity, this parameter should typically be less than one. (optional; default=0.01)}
\item{MaxAtomsP}{Maximum number of atoms in the atomic domain used for the prior of the pattern matrix in the decomposition (see Sibisi and Skilling, 1997).  The default value will typically be sufficient for most applications (optional; default=$2^32$).}
\item{alphaP}{Sparsity parameter reflecting the expected number of atoms per element of the pattern matrix in the decomposition.  To enforce sparsity, this parameter should typically be less than one. (optional; default=0.01)}
\item{SAIter}{Number of burn-in iterations for the MCMC matrix decomposition (optional; default=1000000000)}
\item{iter}{Number of iterations to represent the distribution of amplitude and pattern matrices with the MCMC matrix decomposition (optional; default=500000000)}
\item{thin}{Double whose integer part represents the number of iterations at which the samples are kept and decimal part provides an identifier for the output files from this implementation of GAPS.  If thin is an integer or not specified, this decimal file identifier is assigned randomly.  (optional; default=-1; code assigns number of iterations kept to be iter/10000 and file identifier to be runif(1)) }
\item{verbose}{Boolean which specifies the amount of output to the user about the progress of the program. (optional; default=TRUE)}
\item{keepChain}{Boolean which specifies if chain values of \eqn{{\bf{A}}} and \eqn{{\bf{P}}} are saved in outputDir (optional; default=FALSE).}
}

\details{The decomposition in GAPS is achieved by finding amplitude and pattern matrices (\eqn{{\bf{A}}} and \eqn{{\bf{P}}}, respectively) for which \deqn{{\bf{D}} = {\bf{A}}{\bf{P}} + \Sigma}, where \eqn{\Sigma} is the matrix of uncertainties given by unc.  The matrices \eqn{\bf{A}} and \eqn{\bf{P}} are assumed to have the atomic prior described in Sibisi and Skilling (1997) and are found with MCMC sampling implemented within JAGS.}

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
}

\note{Running GAPS will create the folder ouptutDir, create diagnostic files with \eqn{\chi^{2}} and number of atoms, files with the mean and standard deviation of \eqn{{\bf{A}}} and \eqn{{\bf{P}}}, and optionally values of \eqn{{\bf{A}}} and \eqn{{\bf{P}}} from the MCMC chain.}

\author{Elana J. Fertig \email{ejfertig@jhmi.edu}}

\references{
M. Plummer. (2003) JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling. In K. Hornik, F. Leisch, and A. Zeileis, editors, Proceedings of the Third International Workshop on Distributed Statistical Computing, Vienna, Austria.

S. Sibisi and J. Skilling. (1997) Prior distributions on measure space.  Journal of the Royal Statistical Society, B, 59:217-235.
}

\examples{
\dontrun{
## Load data
data(ModSim)

## Run GAPS matrix decomposition
nIter <- 500000
results <- GAPS(data=ModSim.D, unc=0.01, isPercentError=FALSE,
                numPatterns=3, SAIter=2*nIter, iter = nIter,
                outputDir='ModSimResults')

## Plot the results
plotGAPS(results$Amean, results$Pmean)
}
}

\seealso{\code{\link{CoGAPS}}} 
\keyword{misc}
