context("C++")

test_that("Catch unit tests pass",
{
    data(GIST)

    gistCsvPath <<- system.file("extdata/GIST.csv", package="CoGAPS")
    gistTsvPath <<- system.file("extdata/GIST.tsv", package="CoGAPS")
    gistMtxPath <<- system.file("extdata/GIST.mtx", package="CoGAPS")
    gistGctPath <<- system.file("extdata/GIST.gct", package="CoGAPS")

    run_catch_unit_tests()
    expect(TRUE)
})
