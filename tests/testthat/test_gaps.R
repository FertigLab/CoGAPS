context("GAPS")

test_that("GAPS Simple Simulation", {
    data(SimpSim)
    nIter <- 1000
    results <- gapsRun(SimpSim.D, SimpSim.S, nFactor=3,
            nEquil=nIter, nSample=nIter)
})

test_that("GAPSmap Simple Simulation", {
    data(SimpSim)
    FP <- matrix(SimpSim.P[3, ], nrow=1)
    nIter <- 1000
    results <- gapsMapRun(SimpSim.D, SimpSim.S, FP, nFactor=3,
                          nEquil=nIter, nSample=nIter)
})
