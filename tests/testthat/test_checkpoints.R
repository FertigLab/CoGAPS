context("CoGAPS")

test_that("Checkpoint System",
{
    if (CoGAPS::checkpointsEnabled())
    {
        data(GIST)
        run1 <- CoGAPS(GIST.matrix, checkpointInterval=51,
            checkpointOutFile="test.out", messages=FALSE, nIterations=100)
        run2 <- CoGAPS(GIST.matrix, checkpointInFile="test.out", messages=FALSE,
            nIterations=100)

        print(max(run1@featureLoadings - run2@featureLoadings))

        expect_true(all.equal(run1@featureLoadings, run2@featureLoadings))
        expect_true(all.equal(run1@sampleFactors, run2@sampleFactors))
    }
})