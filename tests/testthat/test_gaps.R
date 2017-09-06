context("GAPS")

test_that("GAPS Simple Simulation", {
    data(SimpSim)
    nIter <- 1000
    expect_failure(expect_error(gapsRun(SimpSim.D, SimpSim.S, nFactor=3,
                                        messages=FALSE)))
})

test_that("GAPSmap Simple Simulation", {
    data(SimpSim)
    FP <- matrix(SimpSim.P[3, ], nrow=1)
    nIter <- 1000
    expect_failure(expect_error(gapsMapRun(SimpSim.D, SimpSim.S, FP,
                                           nFactor=3, messages=FALSE)))
})
