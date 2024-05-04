library('testthat')

test_that("equal A and P dimensions in sparse vs standard", {
    data(GIST)
    res_standard <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100,
                           seed=1, messages=FALSE)
    res_sparse <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100, seed=1,
                         messages=FALSE, sparseOptimization=TRUE)
    expect_equal(dim(res_standard@featureLoadings), dim(res_sparse@featureLoadings))
    expect_equal(dim(res_standard@sampleFactors), dim(res_sparse@sampleFactors))
})

test_that("equal A and P dimensions in sc vs standard", {
    data(GIST)
    res_standard <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100,
                           seed=1, messages=FALSE)
    params <- CogapsParams()
    params <- setDistributedParams(params, nSets=2)
    res_sc <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100, seed=1,
                     messages=FALSE, distributed="single-cell")
    expect_equal(dim(res_standard@featureLoadings), dim(res_sc@featureLoadings))
    expect_equal(dim(res_standard@sampleFactors), dim(res_sc@sampleFactors))
})

test_that("equal A and P dimensions in gw vs standard", {
    data(GIST)
    res_standard <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100,
                           seed=1, messages=FALSE)
    params <- CogapsParams()
    params <- setDistributedParams(params, nSets=2)
    res_gw <- CoGAPS(GIST.data_frame, nPatterns=2, nIterations=100, seed=1,
                     messages=FALSE, distributed="genome-wide")
    expect_equal(dim(res_standard@featureLoadings), dim(res_gw@featureLoadings))
    expect_equal(dim(res_standard@sampleFactors), dim(res_gw@sampleFactors))
})