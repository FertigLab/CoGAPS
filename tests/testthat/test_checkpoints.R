context("CoGAPS")

test_that("Checkpoint System",
{
    data(SimpSim)
    run1 <- CoGAPS(SimpSim.data, checkpointInterval=501,
        checkpointOutFile="test.out", messages=FALSE)
    run2 <- CoGAPS(SimpSim.data, checkpointInFile="test.out", messages=FALSE)

    print(max(run1@featureLoadings - run2@featureLoadings))

    expect_true(all.equal(run1@featureLoadings, run2@featureLoadings))
    expect_true(all.equal(run1@sampleFactors, run2@sampleFactors))
})