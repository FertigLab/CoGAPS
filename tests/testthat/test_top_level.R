context("CoGAPS")

test_that("Valid Top-Level CoGAPS Calls",
{
    data(GIST)
    testDataFrame <- GIST.data_frame
    testMatrix <- GIST.matrix
    #testSummarizedExperiment <- 
    #testSingleCellExperiment

    gistCsvPath <- system.file("extdata/GIST.csv", package="CoGAPS")
    gistTsvPath <- system.file("extdata/GIST.tsv", package="CoGAPS")
    gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")

    # data types
    res <- list()
    res[[1]] <- CoGAPS(testDataFrame, nIterations=100, outputFrequency=50, seed=1)
    res[[2]] <- CoGAPS(testMatrix, nIterations=100, outputFrequency=50, seed=1)
    res[[3]] <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=50, seed=1)
    res[[4]] <- CoGAPS(gistTsvPath, nIterations=100, outputFrequency=50, seed=1)
    res[[5]] <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=50, seed=1)
    
    expect_equal(nrow(res[[1]]@featureLoadings), 1363)
    expect_equal(ncol(res[[1]]@featureLoadings), 7)
    expect_equal(nrow(res[[1]]@sampleFactors), 9)
    expect_equal(ncol(res[[1]]@sampleFactors), 7)
    expect_true(all(sapply(1:4, function(i)
        res[[i]]@featureLoadings == res[[i+1]]@featureLoadings)))
    expect_true(all(sapply(1:4, function(i)
        res[[i]]@sampleFactors == res[[i+1]]@sampleFactors)))

    # transposing data
    res <- list()
    res[[1]] <- CoGAPS(testDataFrame, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1)
    res[[2]] <- CoGAPS(testMatrix, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1)
    res[[3]] <- CoGAPS(gistCsvPath, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1)
    res[[4]] <- CoGAPS(gistTsvPath, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1)
    res[[5]] <- CoGAPS(gistMtxPath, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1)
    
    expect_equal(nrow(res[[1]]@featureLoadings), 9)
    expect_equal(ncol(res[[1]]@featureLoadings), 7)
    expect_equal(nrow(res[[1]]@sampleFactors), 1363)
    expect_equal(ncol(res[[1]]@sampleFactors), 7)
    expect_true(all(sapply(1:4, function(i)
        res[[i]]@featureLoadings == res[[i+1]]@featureLoadings)))
    expect_true(all(sapply(1:4, function(i)
        res[[i]]@sampleFactors == res[[i+1]]@sampleFactors)))

    # passing uncertainty
    #expect_error(CoGAPS(testDataFrame, uncertainty=as.matrix(GIST.S),
    #    nIterations=100, outputFrequency=50, seed=1), NA)    

    # multiple threads
    expect_error(CoGAPS(testDataFrame, nIterations=100, outputFrequency=50,
        seed=1, nThreads=2), NA)
    expect_error(CoGAPS(testDataFrame, nIterations=100, outputFrequency=50,
        seed=1, nThreads=6), NA)
    expect_error(CoGAPS(testDataFrame, nIterations=100, outputFrequency=50,
        seed=1, nThreads=12), NA)

    # distributed CoGAPS 
    expect_error(CoGAPS(testDataFrame, nIterations=100, outputFrequency=50,
        seed=1, distributed="genome-wide"), NA)
    expect_error(CoGAPS(gistTsvPath, nIterations=100, outputFrequency=50,
        seed=1, distributed="genome-wide"), NA)

    expect_error(CoGAPS(testDataFrame, nIterations=100, outputFrequency=50,
        seed=1, distributed="single-cell", singleCell=TRUE,
        transposeData=TRUE), NA)
    expect_error(CoGAPS(gistTsvPath, nIterations=100, outputFrequency=50,
        seed=1, distributed="single-cell", singleCell=TRUE,
        transposeData=TRUE), NA)
})

#test_that("Invalid Top-Level CoGAPS Calls",
#{

#})