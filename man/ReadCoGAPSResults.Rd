\name{ReadCoGAPSResults}
\alias{ReadCoGAPSResults}
\title{Parse results saved to output folder after CoGAPS simulations}

\description{Read all instances of CoGAPS simulations in a specified directory.}

\usage{
  ReadCoGAPSResults(path=getwd(), output.list=TRUE)
}

\arguments{
\item{path}{Directory containing all results from the CoGAPS simulation (optional; defaults to current working directory).}
\item{output.list}{Boolean specifying whether each simulation should be output as a list element or as rows/columns appended to the appropriate matrices (optional; defaults to list output).}
}

\value{
  A list containing:
  \item{A.mean}{Mean amplitude matrix of each simulation as a list (output.list=T) or a matrix (output.list=F).}
  \item{A.sd}{Standard deviation of the amplitude matrix of each simulation as a list (output.list=T) or a matrix (output.list=F).}
  \item{P.mean}{Mean pattern matrix of each simulation as a list (output.list=T) or a matrix (output.list=F).}
  \item{P.sd}{Standard deviation of the pattern matrix of each simulation as a list (output.list=T) or a matrix (output.list=F).}
  \item{M}{Mock data (noise filtered data) for each simulation generated from \code{A.mean \%*\% P.mean} for each simulation as a list (output.list=T) or a matrix (output.list=F).}
}
