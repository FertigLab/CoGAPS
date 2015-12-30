#' CoGAPS: Coordinated Gene Activity in Pattern Sets
#'
#' CoGAPS implements a Bayesian MCMC matrix factorization algorithm,
#' GAPS, and links it to gene set statistic methods to infer
#' biological process activity.  It can be used to perform
#' sparse matrix factorization on any data, and when this
#' data represents biomolecules, to do gene set analysis.
#' \tabular{ll}{
#' Package: \tab CoGAPS\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2014-07-23\cr
#' License: \tab LGPL\cr
#' }
#' @author Maintainer: Elana J. Fertig \email{ejfertig@jhmi.edu},
#'             Michael F. Ochs \email{ochsm@tcnj.edu}
#' @references
#' Fertig EJ, Ding J, Favorov AV, Parmigiani G, Ochs MF.
#' CoGAPS: an R/C++ package to identify patterns and biological
#' process activity in transcriptomic data.
#' Bioinformatics. 2010 Nov 1;26(21):2792-3
#' @docType package
#' @name CoGAPS-package
#' @importFrom Rcpp evalCpp
#' @importFrom gplots heatmap.2 plotCI
#' @importFrom stats variable.names sd update heatmap runif
#' @importFrom graphics matplot title abline close.screen hist legend lines mtext par plot points screen split.screen
#' @importFrom grDevices dev.new dev.off pdf colorRampPalette rainbow
#' @importFrom methods is
#' @importFrom utils read.table write.table
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats cor loess
#' @useDynLib CoGAPS
NULL
