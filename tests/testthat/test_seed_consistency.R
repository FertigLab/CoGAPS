context("CoGAPS")

checkCompare <- function(comp)
{
    if (is(comp, "character"))
    {
        print(comp)
        return(FALSE)
    }
    return(TRUE)    
}

resultsEqual <- function(res1, res2)
{
    checkCompare(all.equal(res1@featureLoadings, res2@featureLoadings, tolerance=0.1)) &
    checkCompare(all.equal(res1@featureStdDev, res2@featureStdDev, tolerance=0.1)) &
    checkCompare(all.equal(res1@sampleFactors, res2@sampleFactors, tolerance=0.1)) &
    checkCompare(all.equal(res1@sampleStdDev, res2@sampleStdDev, tolerance=0.1))
    checkCompare(all.equal(res1@metadata$atomsA, res2@metadata$atomsA)) &
    checkCompare(all.equal(res1@metadata$atomsP, res2@metadata$atomsP))
}

test_that("same seed == same result",
{
    gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")

    # standard cogaps
    res1 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE)
    res2 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE)
    expect_true(resultsEqual(res1, res2))

    # distributed cogaps
    res1 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, distributed="genome-wide")
    res2 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, distributed="genome-wide")
    expect_true(resultsEqual(res1, res2))
    
    # multiple threads, dense sampler
    res1 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, nThreads=1, sparseOptimization=FALSE)
    res2 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, nThreads=3, sparseOptimization=FALSE)
    res3 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, nThreads=6, sparseOptimization=FALSE)
    expect_true(resultsEqual(res1, res2))
    expect_true(resultsEqual(res1, res3))
    expect_true(resultsEqual(res2, res3))

    # multiple threads, sparse sampler
    res1 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, nThreads=1, sparseOptimization=TRUE)
    res2 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, nThreads=3, sparseOptimization=TRUE)
    res3 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, nThreads=6, sparseOptimization=TRUE)
    expect_true(resultsEqual(res1, res2))
    expect_true(resultsEqual(res1, res3))
    expect_true(resultsEqual(res2, res3))
})