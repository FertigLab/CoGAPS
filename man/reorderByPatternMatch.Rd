\name{reorderByPatternMatch}
\alias{reorderByPatternMatch}
\title{Match two sets of patterns found with CoGAPS}

\description{Matches two sets of pattern matrices (of the same size) found with CoGAPS.  Matches are identified by finding identifying subsequently decreasing correlations between patterns in the respective matrices.}

\usage{
  reorderByPatternMatch(P, matchTo) 
}

\arguments{
\item{P}{Pattern matrix for which rows will be arranged to match the matrix in matchTo}
\item{matchTo}{Pattern matrix to which P is matched.}
}

\value{Pattern matrix derived from reordering columns of P}

