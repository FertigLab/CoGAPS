context("CoGAPS")

test_that("Npatterns are required input",
{
    expect_error(CoGAPS(nIterations=100, messages=FALSE))
    expect_no_error(CoGAPS(nPatterns=7, nIterations=100, messages=FALSE))
})