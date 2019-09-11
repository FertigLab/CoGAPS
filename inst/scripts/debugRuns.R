library(CoGAPS)

# load and create all data used for the runs
data(GIST)
GISTPathCsv <- system.file("extdata/GIST.csv", package="CoGAPS")
GISTPathTsv <- system.file("extdata/GIST.csv", package="CoGAPS")
GISTPathGct <- system.file("extdata/GIST.csv", package="CoGAPS")
GISTPathMtx <- system.file("extdata/GIST.csv", package="CoGAPS")

# default parameters for all data types
CoGAPS(GIST.matrix, nIterations=1000)
CoGAPS(GIST.data_frame, nIterations=1000)
CoGAPS(GISTPathCsv, nIterations=1000)
CoGAPS(GISTPathTsv, nIterations=1000)
CoGAPS(GISTPathGct, nIterations=1000)
CoGAPS(GISTPathMtx, nIterations=1000)

# all sampler types
CoGAPS(GIST.matrix, sparseOptimization=FALSE, asynchronousUpdates=FALSE, nIterations=1000)
CoGAPS(GIST.matrix, sparseOptimization=TRUE, asynchronousUpdates=FALSE, nIterations=1000)
CoGAPS(GIST.matrix, sparseOptimization=FALSE, asynchronousUpdates=TRUE, nIterations=1000)
CoGAPS(GIST.matrix, sparseOptimization=TRUE, asynchronousUpdates=TRUE, nIterations=1000)

# all data types and file readers
CoGAPS(GISTPathCsv, sparseOptimization=FALSE, nIterations=1000)
CoGAPS(GISTPathCsv, sparseOptimization=TRUE, nIterations=1000)
CoGAPS(GISTPathMtx, sparseOptimization=FALSE, nIterations=1000)
CoGAPS(GISTPathMtx, sparseOptimization=TRUE, nIterations=1000)

# multiple threads for all data types
CoGAPS(GISTPathCsv, sparseOptimization=FALSE, nThreads=2, nIterations=1000)
CoGAPS(GISTPathCsv, sparseOptimization=TRUE, nThreads=2, nIterations=1000)

# distributed version
params <- CogapsParams()
params <- setDistributedParams(params, nSets=4)
CoGAPS(GISTPathCsv, params, distributed="genome-wide", nIterations=1000)
CoGAPS(GISTPathCsv, params, distributed="single-cell", transposeData=TRUE, singleCell=TRUE, nIterations=1000)
