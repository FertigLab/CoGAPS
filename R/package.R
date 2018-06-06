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
#' Version: \tab 2.99.0\cr
#' Date: \tab 2018-01-24\cr
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
#' @importFrom utils read.table write.table str
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats cor loess as.dist cutree as.hclust complete.cases D hclust weighted.mean
#' @import foreach
#' @importFrom cluster agnes 
#' @import doParallel
#' @import shiny
#' @importFrom reshape2 melt
#' @import ggplot2
#' @useDynLib CoGAPS
NULL

#' Simulated data
#' @docType data
#' @name SimpSim.D
#' @usage SimpSim.D
NULL

#' Simulated data
#' @docType data
#' @name SimpSim.S
#' @usage SimpSim.S
NULL

#' Simulated data
#' @docType data
#' @name SimpSim.A
#' @usage SimpSim.A
NULL

#' Simulated data
#' @docType data
#' @name SimpSim.P
#' @usage SimpSim.P
NULL

#' CoGAPS result from running on simulated data
#' @docType data
#' @name SimpSim.result
#' @usage SimpSim.result
NULL

#' Sample GIST gene expression data from Ochs et al. (2009)
#' @docType data
#' @name GIST.D
#' @usage GIST.D
NULL

#' Sample GIST gene expression data from Ochs et al. (2009)
#' @docType data
#' @name GIST.S
#' @usage GIST.S
NULL

#' Simulated dataset to quantify gene set membership in the GIST dataset 
#' @docType data
#' @name GIST.GeneSets
#' @usage GIST.GeneSets
NULL

#' CoGAPS result from running on GIST dataset
#' @docType data
#' @name GIST.result
#' @usage GIST.result
NULL

#' Simulated dataset to quantify gene set membership.
#' @docType data
#' @name GSets
#' @usage GSets
NULL

#' Gene sets defined by transcription factors defined from TRANSFAC.
#' @docType data
#' @name tf2ugFC
#' @usage tf2ugFC
NULL