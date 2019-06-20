context("CoGAPS")

resultsEqual <- function(res1, res2)
{
    all(res1@featureLoadings == res2@featureLoadings) &
    all(res1@featureStdDev == res2@featureStdDev) &
    all(res1@sampleFactors == res2@sampleFactors) &
    all(res1@sampleStdDev== res2@sampleStdDev)
}

test_that("same seed == same result",
{
    gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")

    # standard cogaps
    res1 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE)
    res2 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE)
    expect_true(resultsEqual(res1, res2))

    # distributed cogaps
    res1 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, distributed="genome-wide")
    res2 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, distributed="genome-wide")
    expect_true(resultsEqual(res1, res2))
    
    # multiple threads
    res1 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nThreads=1)
    res2 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nThreads=3)
    res3 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=100, seed=42,
        messages=FALSE, nThreads=6)
    expect_true(resultsEqual(res1, res2))
    expect_true(resultsEqual(res1, res3))
    expect_true(resultsEqual(res2, res3))
})