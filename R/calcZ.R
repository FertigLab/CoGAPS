# calcZ: function to generate Z-score matrix
# History: v1.0 EJF original CoGAPS

# Inputs: meanMat - matrix of mean values
#         sdMat - matrix of standard deviations

# Output: matrix of Z-scores

#'\code{calcZ} calculates the Z-score for each element based
#'on input mean and standard deviation matrices
#'
#'@param meanMat matrix of mean values
#'@param sdMat matrix of standard deviation values
#'@export

calcZ <- function (meanMat, sdMat) {

  # compute the z score for each gene's association to each pattern

  # find matrix dimensions
  nrows <- dim(meanMat)[1]
  ncols <- dim(meanMat)[2]

  check <- dim(sdMat)[1]
  if (nrows != check) {
    stop("Number of rows in the mean and standard deviation of A do not agree.")
  }

  check <- dim(sdMat)[2]
  if (ncols != check) {
    stop("Number of columns in the mean and standard deviation of A do not agree.")
  }

  # compute the matrix of z scores
  zMat <- meanMat/sdMat
  rownames(zMat) <- rownames(meanMat)
  colnames(zMat) <- colnames(meanMat)

  return(zMat)
}
