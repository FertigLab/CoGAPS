context("Snapshots")

test_that("Distinct Snapshots", {
    data(SimpSim)
    x <- gapsRun(SimpSim.D,SimpSim.S,nFactor = 3)
    expect_false(identical(x$ASnapshots[, , 1], x$ASnapshots[, , 100]))
    expect_false(identical(x$PSnapshots[, , 1], x$PSnapshots[, , 100]))
})
