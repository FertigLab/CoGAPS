context("C++")

test_that("Catch unit tests pass",
{
    data(SimpSim)
    data(GIST_TS_20084)
    run_catch_unit_tests()
})
