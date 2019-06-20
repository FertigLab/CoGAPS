library(CoGAPS)

# load and create all data used for the runs
data(GIST)
GISTPathCsv <- system.file("extdata/GIST.csv", package="CoGAPS")
GISTPathTsv <- system.file("extdata/GIST.csv", package="CoGAPS")
GISTPathGct <- system.file("extdata/GIST.csv", package="CoGAPS")
GISTPathMtx <- system.file("extdata/GIST.csv", package="CoGAPS")

# default parameters for all data types
CoGAPS(GIST.matrix)
CoGAPS(GIST.data_frame)
CoGAPS(GISTPathCsv)
CoGAPS(GISTPathTsv)
CoGAPS(GISTPathGct)
CoGAPS(GISTPathMtx)

# all sampler types
CoGAPS(GIST.matrix, sparseOptimization=FALSE, asynchronousUpdates=FALSE)
CoGAPS(GIST.matrix, sparseOptimization=TRUE, asynchronousUpdates=FALSE)
CoGAPS(GIST.matrix, sparseOptimization=FALSE, asynchronousUpdates=TRUE)
CoGAPS(GIST.matrix, sparseOptimization=TRUE, asynchronousUpdates=TRUE)

# all data types and file readers
CoGAPS(GISTPathCsv, sparseOptimization=FALSE)
CoGAPS(GISTPathCsv, sparseOptimization=TRUE)
CoGAPS(GISTPathMtx, sparseOptimization=FALSE)
CoGAPS(GISTPathMtx, sparseOptimization=TRUE)


