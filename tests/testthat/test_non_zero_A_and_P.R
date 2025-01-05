context("CoGAPS")


test_that("CoGAPS returns non-zero-populated patterns and loadings",
{
    data(GIST)
    testDataFrame <- GIST.data_frame
    nPat<-5
    res <- CoGAPS(testDataFrame, nPatterns=nPat, nIterations=100, outputFrequency=50, seed=42, messages=TRUE)
    expect_equal(nrow(res@featureLoadings), 1363)
    expect_equal(ncol(res@featureLoadings), nPat)
    expect_equal(nrow(res@sampleFactors), 9)
    expect_equal(ncol(res@sampleFactors), nPat)
    expect_true(all(apply(res@sampleFactors,2,sum) != 0))
    expect_true(all(apply(res@featureLoadings,2,sum) != 0))
})

