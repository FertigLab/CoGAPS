context("Snapshots")

test_that("Distinct gaps Snapshots", {
    data(SimpSim)
    x <- gapsRun(SimpSim.D, SimpSim.S, nFactor=3, messages=FALSE)
    expect_false(identical(x$ASnapshots[, , 1], x$ASnapshots[, , 100]))
    expect_false(identical(x$PSnapshots[, , 1], x$PSnapshots[, , 100]))
})

test_that("Distinct gapsMap Snapshots", {
    data(SimpSim)
    FP <- matrix(SimpSim.P[3, ], nrow=1)
    x <- gapsMapRun(SimpSim.D, SimpSim.S, FP, nFactor=3, messages=FALSE)
    expect_false(identical(x$ASnapshots[, , 1], x$ASnapshots[, , 100]))
    expect_false(identical(x$PSnapshots[, , 1], x$PSnapshots[, , 100]))
})
