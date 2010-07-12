\name{plotGAPS}
\alias{plotGAPS}
\title{Plotter for GAPS decomposition results}

\description{Plots the A and P matrices obtained from the GAPS matrix decomposition.}

\usage{
  plotGAPS(A, P, outputPDF="")
}

\arguments{
\item{A}{The amplitude matrix \eqn{{\bf{A}}} obtained from GAPS.} 
\item{P}{The pattern matrix \eqn{{\bf{P}}} obtained from GAPS.}
\item{outputPDF}{Name of an pdf file to which the results will be output.  (Optional; default="" will output plots to screen).}
}

\note{
If the plot option is true in \code{\link{CoGAPS}}, this function will be called automatically to plot results to the screen.
}

\author{Elana J. Fertig \email{efertig@jhmi.edu}}

\seealso{\code{\link{CoGAPS}}}
