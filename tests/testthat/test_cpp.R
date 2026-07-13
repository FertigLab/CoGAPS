context("C++")

test_that("Catch unit tests pass",
{
    # the C++ file-parser tests read these data paths from the global environment
    gistCsvPath <<- system.file("extdata/GIST.csv", package="CoGAPS")
    gistTsvPath <<- system.file("extdata/GIST.tsv", package="CoGAPS")
    gistMtxPath <<- system.file("extdata/GIST.mtx", package="CoGAPS")
    gistGctPath <<- system.file("extdata/GIST.gct", package="CoGAPS")

    # run_catch_unit_tests() runs the Catch2 C++ unit test suite and returns the
    # number of failed assertions (0 == all pass). It is an internal function, so
    # it must be reached through the ::: operator. Asserting == 0 is what actually
    # propagates a C++ test failure up to `R CMD check`.
    expect_equal(CoGAPS:::run_catch_unit_tests(), 0L)
})
