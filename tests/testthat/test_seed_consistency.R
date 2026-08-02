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
    checkCompare(all.equal(res1@loadingStdDev, res2@loadingStdDev, tolerance=0.1)) &
    checkCompare(all.equal(res1@sampleFactors, res2@sampleFactors, tolerance=0.1)) &
    checkCompare(all.equal(res1@factorStdDev, res2@factorStdDev, tolerance=0.1))
    checkCompare(all.equal(res1@metadata$atomsA, res2@metadata$atomsA)) &
    checkCompare(all.equal(res1@metadata$atomsP, res2@metadata$atomsP))
}

test_that("same seed == same result",
{
    gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")

    # standard cogaps
    res1 <- CoGAPS(gistMtxPath, nPatterns=7, nIterations=100, outputFrequency=10,
                   seed=42, messages=FALSE)
    res2 <- CoGAPS(gistMtxPath, nPatterns=7, nIterations=100, outputFrequency=10,
                   seed=42, messages=FALSE)
    expect_true(resultsEqual(res1, res2))

    # distributed cogaps
    res1 <- CoGAPS(gistMtxPath, nPatterns=7, nIterations=100, outputFrequency=10,
                   seed=42, messages=FALSE, distributed="genome-wide")
    res2 <- CoGAPS(gistMtxPath, nIterations=100, outputFrequency=10, seed=42,
                   nPatterns=7, messages=FALSE, distributed="genome-wide")
    expect_true(resultsEqual(res1, res2))
    
    # seed consistency, dense sampler
    res1 <- CoGAPS(gistMtxPath, nPatterns=7, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, sparseOptimization=FALSE)
    res2 <- CoGAPS(gistMtxPath, nPatterns=7, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, sparseOptimization=FALSE)
    expect_true(resultsEqual(res1, res2))

    # seed consistency, sparse sampler
    res1 <- CoGAPS(gistMtxPath, nPatterns=7, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, sparseOptimization=TRUE)
    res2 <- CoGAPS(gistMtxPath, nPatterns=7, nIterations=100, outputFrequency=10, seed=42,
        messages=FALSE, sparseOptimization=TRUE)
    expect_true(resultsEqual(res1, res2))
})