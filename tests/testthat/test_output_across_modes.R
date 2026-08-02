library('testthat')

test_that("equal A and P dimensions in sparse vs standard", {
    data(GIST)
    res_standard <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100,
                           seed=1, messages=FALSE)
    res_sparse <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100, seed=1,
                         messages=FALSE, sparseOptimization=TRUE)
    expect_equal(dim(res_standard@featureLoadings),
                 dim(res_sparse@featureLoadings))
    expect_equal(dim(res_standard@sampleFactors), 
                 dim(res_sparse@sampleFactors))
})

# distributed CoGAPS needs on-disk data (mtx/tsv/csv/gct); passing an in-memory
# matrix warns. Use the packaged GIST.mtx, which is the same data as GIST.data_frame.
test_that("equal A and P dimensions in sc vs standard", {
    data(GIST)
    gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")
    res_standard <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100,
                           seed=1, messages=FALSE)
    res_sc <- CoGAPS(gistMtxPath, nPatterns=2, nIterations=100, seed=1,
                     messages=FALSE, distributed="single-cell")
    expect_equal(dim(res_standard@featureLoadings), dim(res_sc@featureLoadings))
    expect_equal(dim(res_standard@sampleFactors), dim(res_sc@sampleFactors))
})

test_that("equal A and P dimensions in gw vs standard", {
    data(GIST)
    gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")
    res_standard <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100,
                           seed=1, messages=FALSE)
    res_gw <- CoGAPS(gistMtxPath, nPatterns=2, nIterations=100, seed=1,
                     messages=FALSE, distributed="genome-wide")
    expect_equal(dim(res_standard@featureLoadings), dim(res_gw@featureLoadings))
    expect_equal(dim(res_standard@sampleFactors), dim(res_gw@sampleFactors))
})