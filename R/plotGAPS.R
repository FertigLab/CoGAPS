#' Plot Decomposed A and P Matrices
#'
#' @details plots the output A and P matrices as a
#' heatmap and line plot respectively
#' @param A the mean A matrix
#' @param P the mean P matrix
#' @param outputPDF optional root name for PDF output, if
#' not specified, output goes to screen
#' @return plot
#' @examples
#' data(SimpSim)
#' plotGAPS(SimpSim.result$Amean, SimpSim.result$Pmean)
#' @export
plotGAPS <- function(A, P)
{
    plotP(P)
    heatmap(A, Rowv=NA, Colv=NA)
}
