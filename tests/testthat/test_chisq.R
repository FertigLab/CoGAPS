test_that('chi-square reported by CoGAPS mathches manually calculated (w/uncertainty)',{
    data(GIST)
    data <- GIST.data_frame
    unc <- 0.1*as.matrix(data)
    res <- CoGAPS(data, nIterations=1000, uncertainty=unc, nPatterns=3,
                  seed=1, messages=FALSE, sparseOptimization=FALSE)
    reported <- getMeanChiSq(res)

    A <- getAmplitudeMatrix(res)
    P <- getPatternMatrix(res)
    M <- A %*% t(P)

    calculated <- sum(((data - M)/unc)^2)
    #Set the tolerance to float precision times matrix size
    sprintf("difference: %.17f", abs(reported - calculated))
    expect_equal(reported, calculated, tolerance = (1e-7)*prod(dim(data)))
})
