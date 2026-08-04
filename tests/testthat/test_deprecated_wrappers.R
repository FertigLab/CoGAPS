context("CoGAPS")

# scCoGAPS() and GWCoGAPS() are deprecated wrappers around
# CoGAPS(..., distributed=). They are still exported, so they have to keep
# working. Neither was called anywhere in the test suite or in a runnable
# example, which is how they came to be broken unnoticed: once nPatterns became
# mandatory, their `params=new("CogapsParams")` default could no longer be
# evaluated, and every call failed -- including calls that passed nPatterns,
# because the wrapper touches `params` before forwarding `...`.

# nSets cannot be passed through `...` -- it has to go through
# setDistributedParams() -- so the no-params cases rely on the default.
gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")

test_that("scCoGAPS accepts nPatterns without an explicit params object",
{
    expect_warning(res <- CoGAPS::scCoGAPS(gistMtxPath, nPatterns=2,
        nIterations=50, seed=1, messages=FALSE,
        BPPARAM=BiocParallel::SerialParam()), "deprecated")
    expect_true(is(res, "CogapsResult"))
    expect_equal(ncol(res@featureLoadings), 2)
})

test_that("GWCoGAPS accepts nPatterns without an explicit params object",
{
    expect_warning(res <- CoGAPS::GWCoGAPS(gistMtxPath, nPatterns=2,
        nIterations=50, seed=1, messages=FALSE,
        BPPARAM=BiocParallel::SerialParam()), "deprecated")
    expect_true(is(res, "CogapsResult"))
    expect_equal(ncol(res@featureLoadings), 2)
})

test_that("the wrappers still accept an explicit params object",
{
    params <- CogapsParams(nPatterns=2)
    params <- setDistributedParams(params, nSets=2)

    expect_warning(res <- CoGAPS::scCoGAPS(gistMtxPath, params=params,
        nIterations=50, seed=1, messages=FALSE,
        BPPARAM=BiocParallel::SerialParam()), "deprecated")
    expect_equal(ncol(res@featureLoadings), 2)
})

test_that("the wrappers set the distributed mode they are named after",
{
    params <- CogapsParams(nPatterns=2)
    params <- setDistributedParams(params, nSets=2)

    expect_warning(sc <- CoGAPS::scCoGAPS(gistMtxPath, params=params,
        nIterations=50, seed=1, messages=FALSE,
        BPPARAM=BiocParallel::SerialParam()), "deprecated")
    expect_equal(sc@metadata$params@distributed, "single-cell")

    expect_warning(gw <- CoGAPS::GWCoGAPS(gistMtxPath, params=params,
        nIterations=50, seed=1, messages=FALSE,
        BPPARAM=BiocParallel::SerialParam()), "deprecated")
    expect_equal(gw@metadata$params@distributed, "genome-wide")
})
