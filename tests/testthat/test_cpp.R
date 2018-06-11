context("C++")

test_that("Catch unit tests pass",
{
    data(SimpSim)
    data(GIST)

    run_catch_unit_tests()
})
