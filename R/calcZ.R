#' Compute Z-Score Matrix
#'
#' @details calculates the Z-score for each element based on input mean
#' and standard deviation matrices
#' @param meanMat matrix of mean values
#' @param sdMat matrix of standard deviation values
#' @return matrix of z-scores
#' @examples
#' data(SimpSim)
#' calcZ(SimpSim.result$Amean, SimpSim.result$Asd)
#' @export
calcZ <- function(meanMat, sdMat)
{
    if (nrow(meanMat) != nrow(sdMat) | ncol(meanMat) != ncol(sdMat))
        stop("meanMat and sdMat dimensions don't match")

    zMat <- meanMat / sdMat
    rownames(zMat) <- rownames(meanMat)
    colnames(zMat) <- colnames(meanMat)
    return(zMat)
}
