context("Seed")

test_that("Same seed, same output with GAPS Simple Simulation", {
    # load data
    data(SimpSim)

    # number of burn-ins and posterior samples
    nIter <- 5000

    # two runs with same seed
    results1 <- gapsRun(SimpSim.D, SimpSim.S, nFactor=3,
                        nEquil=nIter, nSample=nIter, seed=42)
    results2 <- gapsRun(SimpSim.D, SimpSim.S, nFactor=3,
                        nEquil=nIter, nSample=nIter, seed=42)

    # one run with different seed
    results3 <- gapsRun(SimpSim.D, SimpSim.S, nFactor=3,
                        nEquil=nIter, nSample=nIter, seed=8675309)

    # check that same seed results in identical outputs
    expect_identical(results1, results2)

    # check that different seeds yield different outputs
    expect_false(identical(results1, results3))
})
