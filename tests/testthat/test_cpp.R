context("C++")

test_that("Catch unit tests pass",
{
    data(SimpSim)
    data(GIST)

    gistMtxPath <- system.file("data/GIST.mtx", package="CoGAPS")

    run_catch_unit_tests()
})
