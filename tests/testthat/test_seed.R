context("Seed")

test_that("Same seed, same output with GAPS Simple Simulation", {
    data(SimpSim)
    nIter <- 5000
    results <- gapsRun(SimpSim.D, SimpSim.S, nFactor=3,
            nEquil=nIter, nSample=nIter, seed=42)
})
