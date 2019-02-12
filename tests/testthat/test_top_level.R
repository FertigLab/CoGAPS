context("CoGAPS")

no_na_in_result <- function(gapsResult)
{
    sum(is.na(gapsResult@featureLoadings)) +
        sum(is.na(gapsResult@featureStdDev)) +
        sum(is.na(gapsResult@sampleFactors)) +
        sum(is.na(gapsResult@sampleStdDev)) == 0
}

test_that("Valid Top-Level CoGAPS Calls",
{
    data(GIST)
    testDataFrame <- GIST.data_frame
    testMatrix <- GIST.matrix

    gistCsvPath <- system.file("extdata/GIST.csv", package="CoGAPS")
    gistTsvPath <- system.file("extdata/GIST.tsv", package="CoGAPS")
    gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")
    gistGctPath <- system.file("extdata/GIST.gct", package="CoGAPS")

    # data types
    res <- list()
    res[[1]] <- CoGAPS(testDataFrame, nIterations=100, outputFrequency=50, seed=1, messages=FALSE)
    res[[2]] <- CoGAPS(testMatrix, nIterations=100, outputFrequency=50, seed=1, messages=FALSE)
    res[[3]] <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=50, seed=1, messages=FALSE)
    res[[4]] <- CoGAPS(gistTsvPath, nIterations=100, outputFrequency=50, seed=1, messages=FALSE)
    res[[5]] <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=50, seed=1, messages=FALSE)
    res[[6]] <- CoGAPS(gistGctPath, nIterations=100, outputFrequency=50, seed=1, messages=FALSE)
    expect_true(all(sapply(res, no_na_in_result)))
    
    expect_equal(nrow(res[[1]]@featureLoadings), 1363)
    expect_equal(ncol(res[[1]]@featureLoadings), 7)
    expect_equal(nrow(res[[1]]@sampleFactors), 9)
    expect_equal(ncol(res[[1]]@sampleFactors), 7)
    expect_true(all(sapply(1:5, function(i)
        res[[i]]@featureLoadings == res[[i+1]]@featureLoadings)))
    expect_true(all(sapply(1:5, function(i)
        res[[i]]@sampleFactors == res[[i+1]]@sampleFactors)))

    # transposing data
    res <- list()
    res[[1]] <- CoGAPS(testDataFrame, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE)
    res[[2]] <- CoGAPS(testMatrix, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE)
    res[[3]] <- CoGAPS(gistCsvPath, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE)
    res[[4]] <- CoGAPS(gistTsvPath, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE)
    res[[5]] <- CoGAPS(gistMtxPath, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE)
    res[[6]] <- CoGAPS(gistGctPath, transposeData=TRUE, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE)
    expect_true(all(sapply(res, no_na_in_result)))
    
    expect_equal(nrow(res[[1]]@featureLoadings), 9)
    expect_equal(ncol(res[[1]]@featureLoadings), 7)
    expect_equal(nrow(res[[1]]@sampleFactors), 1363)
    expect_equal(ncol(res[[1]]@sampleFactors), 7)
    expect_true(all(sapply(1:5, function(i)
        res[[i]]@featureLoadings == res[[i+1]]@featureLoadings)))
    expect_true(all(sapply(1:5, function(i)
        res[[i]]@sampleFactors == res[[i+1]]@sampleFactors)))

    # passing uncertainty
    expect_error(res <- CoGAPS(testDataFrame, uncertainty=as.matrix(GIST.uncertainty),
        nIterations=100, outputFrequency=50, seed=1, messages=FALSE), NA)    
    expect_true(no_na_in_result(res))

    # multiple threads
    expect_error(res <- CoGAPS(testDataFrame, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE, nThreads=2), NA)
    expect_true(no_na_in_result(res))

    expect_error(res <- CoGAPS(testDataFrame, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE, nThreads=6), NA)
    expect_true(no_na_in_result(res))

    expect_error(res <- CoGAPS(testDataFrame, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE, nThreads=12), NA)
    expect_true(no_na_in_result(res))

    # genome-wide CoGAPS 
    expect_error(res <- CoGAPS(testDataFrame, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE, distributed="genome-wide"), NA)
    expect_true(no_na_in_result(res))

    expect_equal(nrow(res@featureLoadings), 1363)
    expect_equal(nrow(res@sampleFactors), 9)

    expect_error(res <- CoGAPS(gistTsvPath, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE, distributed="genome-wide"), NA)
    expect_true(no_na_in_result(res))

    expect_equal(nrow(res@featureLoadings), 1363)
    expect_equal(nrow(res@sampleFactors), 9)

    # single-cell CoGAPS
    expect_error(res <- CoGAPS(testDataFrame, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE, distributed="single-cell", singleCell=TRUE,
        transposeData=TRUE), NA)
    expect_true(no_na_in_result(res))

    expect_equal(nrow(res@featureLoadings), 9)
    expect_equal(nrow(res@sampleFactors), 1363)

    expect_error(res <- CoGAPS(gistTsvPath, nIterations=100,
        outputFrequency=50, seed=1, messages=FALSE, distributed="single-cell", singleCell=TRUE,
        transposeData=TRUE), NA)
    expect_true(no_na_in_result(res))

    expect_equal(nrow(res@featureLoadings), 9)
    expect_equal(nrow(res@sampleFactors), 1363)

    # test same seed == same result
    res1 <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE)
    res2 <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE)
    expect_true(all(res1@featureLoadings == res2@featureLoadings))
    expect_true(all(res1@featureStdDev == res2@featureStdDev))
    expect_true(all(res1@sampleFactors == res2@sampleFactors))
    expect_true(all(res1@sampleStdDev== res2@sampleStdDev))

    # test same seed == same result for distributed
    res1 <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, distributed="genome-wide")
    res2 <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, distributed="genome-wide")
    expect_true(all(res1@featureLoadings == res2@featureLoadings))
    expect_true(all(res1@featureStdDev == res2@featureStdDev))
    expect_true(all(res1@sampleFactors == res2@sampleFactors))
    expect_true(all(res1@sampleStdDev== res2@sampleStdDev))

    # test same seed == same result for async
    res1 <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nThreads=1)
    res2 <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nThreads=4)
    expect_true(all(res1@featureLoadings == res2@featureLoadings))
    expect_true(all(res1@featureStdDev == res2@featureStdDev))
    expect_true(all(res1@sampleFactors == res2@sampleFactors))
    expect_true(all(res1@sampleStdDev== res2@sampleStdDev))

    # test same seed == same result for different number of threads
    res1 <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nThreads=1)
    res2 <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nThreads=3)
    res3 <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nThreads=6)

    expect_true(all(res1@featureLoadings == res2@featureLoadings))
    expect_true(all(res1@featureStdDev == res2@featureStdDev))
    expect_true(all(res1@sampleFactors == res2@sampleFactors))
    expect_true(all(res1@sampleStdDev== res2@sampleStdDev))

    expect_true(all(res3@featureLoadings == res2@featureLoadings))
    expect_true(all(res3@featureStdDev == res2@featureStdDev))
    expect_true(all(res3@sampleFactors == res2@sampleFactors))
    expect_true(all(res3@sampleStdDev== res2@sampleStdDev))

    # test running with fixed matrix
    nPat <- 3
    fixedA <- matrix(runif(nrow(testMatrix) * nPat, 1, 10), ncol=nPat)
    fixedP <- matrix(runif(ncol(testMatrix) * nPat, 1, 10), ncol=nPat)
    res <- CoGAPS(testMatrix, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nPatterns=nPat, whichMatrixFixed="A",
        fixedPatterns=fixedA)
    
    expect_true(all(dim(res@featureLoadings) == dim(fixedA)))
    for (i in 1:ncol(fixedA))
        fixedA[,i] <- fixedA[,i] * (res@featureLoadings[1,i] / fixedA[1,i])
    all.equal(unname(res@featureLoadings), fixedA, tolerance=0.001)

    res <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nPatterns=nPat, whichMatrixFixed="P",
        fixedPatterns=fixedP)

    expect_true(all(dim(res@sampleFactors) == dim(fixedP)))
    for (i in 1:ncol(fixedP))
        fixedP[,i] <- fixedP[,i] * (res@sampleFactors[1,i] / fixedP[1,i])
    all.equal(unname(res@sampleFactors), fixedP, tolerance=0.001)


    # make sure that "none" gets converted to NULL for distributed
    res <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nPatterns=3, distributed="none")
    expect_true(is.null(res@metadata$params@distributed))    

    params <- new("CogapsParams")
    params <- setParam(params, "distributed", "none")
    res <- CoGAPS(gistCsvPath, params, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nPatterns=3)
    expect_true(is.null(res@metadata$params@distributed))

    # test subsetting indices
    res <- CoGAPS(gistCsvPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nPatterns=3, subsetIndices=1:100, subsetDim=1)
    expect_true(nrow(res@featureLoadings) == 100)
})

