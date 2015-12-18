\name{plotSmoothPatterns}
\alias{plotSmoothPatterns}
\title{Plot loess smoothed CoGAPS patterns}

\description{Plots the sampled mean value of the pattern matrix \eqn{{\bf{P}}} obtained from the CoGAPS matrix factorization vs. a specificed X value for each sample in the columns of \eqn{{\bf{P}}}. Lines plot loess normalized values of \eqn{{\bf{P}}} vs specified X variables.}

\usage{
  plotSmoothPatterns(P, x=NULL, breaks=NULL, breakStyle=T, orderP=!all(is.null(x)), plotPTS=F, pointCol='black', lineCol='grey', add=F, ...)
}

\arguments{
\item{P}{A [p, M] pattern matrix (P.mean) obtained from the CoGAPS matrix factorization.}
\item{x}{A [M, 1] matrix of values for the X axis for each of the corresponding M columns of P.  (Optional: Default: x=1:M)}
\item{breaks}{A vector of X values at which breaks in plotting should occur. Loess lines fit to data will start and stop at breaks. (Optional: Default: no breaks).  May also be specified as an integer to determine the number of equal groups into which to divide the data.}
\item{breakStyle}{A logical vector. If TRUE, the corresponding break will start a new plot on the row for each pattern. If FALSE, a vertical line will demarcate the break point. (Optional: Defaults to all hard breaks). Note, if one logical value is used, that value will determine the break type at each break point.}
\item{orderP}{A logical value.  If TRUE, vertical ordering of patterns will be determined in order of the value of x at which they peak.  If FALSE, vertical ordering will be determined by the rows in the P matrix. (Optional: Default: FALSE)}
\item{plotPTS}{A logical value. If TRUE, plot will include points for each value of the P matrix in addition to the loess smoothed curve.  If FALSE, only the loess smoothed values of P will be plotted. (Optional: Default: FALSE)}
\item{pointCol}{Color of points of the P matrix plotted when plotPTS=TRUE. (Optional: Default: black)}
\item{lineCol}{Color of loess smoothed values of the P matrix. (Optional: Default: grey)}
\item{add}{A logical value.  If TRUE, plot will be added to existing graphics device.  If FALSE, will create a new graphics device. (Optional: Default: FALSE)}
\item{...}{Additional arguments to plotting functions.}
}

\author{Genevieve Stein-O'Brien \email{gsteino1@jhmi.edu}}

\seealso{\code{\link{CoGAPS}}}

\examples{ 
\dontrun{
# create simulated data
P <- rbind(1:10 + rnorm(10), seq(from=10,to=1) + rnorm(10))

# saved as PDF since figure margins are often too large for the null device with this function
# and the null device may also have trouble with the overlay
pdf('Test.pdf', width=10)
plotSmoothPatterns(P=P, x=rep(seq(from=1,to=10,by=2),each=2), breaks=3, breakStyle=c(F,T,T), plotPTS=T)

# demonstrating the overlay of the plot
plotSmoothPatterns(P=P, x=rep(seq(from=1,to=10,by=2),each=2), breaks=c(0.992,3.660,6.340,9.010),   breakStyle=c(F,T,T), plotPTS=F, add=T, lineCol='red')
dev.off()
}
}

