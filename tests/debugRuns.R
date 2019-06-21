library(CoGAPS)

# load and create all data used for the runs
data(GIST)
GISTPathCsv <- system.file("extdata/GIST.csv", package="CoGAPS")
GISTPathTsv <- system.file("extdata/GIST.csv", package="CoGAPS")
GISTPathGct <- system.file("extdata/GIST.csv", package="CoGAPS")
GISTPathMtx <- system.file("extdata/GIST.csv", package="CoGAPS")

# failing tests
CoGAPS(GISTPathMtx, nIterations=1000, outputFrequency=1000, seed=42,
    messages=TRUE, nThreads=4, sparseOptimization=FALSE)
CoGAPS(GISTPathMtx, nIterations=1000, outputFrequency=1000, seed=42,
    messages=TRUE, nThreads=4, sparseOptimization=FALSE)
CoGAPS(GISTPathMtx, nIterations=1000, outputFrequency=1000, seed=42,
    messages=TRUE, nThreads=1, sparseOptimization=TRUE)
CoGAPS(GISTPathMtx, nIterations=1000, outputFrequency=1000, seed=42,
    messages=TRUE, nThreads=1, sparseOptimization=TRUE)
stop()

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

# multiple threads for all data types
CoGAPS(GISTPathCsv, sparseOptimization=FALSE, nThreads=2)
CoGAPS(GISTPathCsv, sparseOptimization=TRUE, nThreads=2)



