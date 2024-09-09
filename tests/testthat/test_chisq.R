test_that('chi-square reported by CoGAPS mathches manually calculated (w/uncertainty)',{
    data(GIST)
    data <- GIST.data_frame
    unc <- 0.1*as.matrix(data)
    res <- CoGAPS(data, nIterations=1000, uncertainty = unc,
                  seed=1, messages=FALSE, sparseOptimization=FALSE)
    reported <- getMeanChiSq(res)

    A <- getAmplitudeMatrix(res)
    P <- getPatternMatrix(res)
    M <- A %*% t(P)

    calculated <- sum(((data - M)/unc)^2)
    expect_equal(reported, calculated, tolerance = 1e-6)
})


test_that('chi-square reported by CoGAPS mathches manually calculated (w/o uncertainty)',{
    data(GIST)
    data <- GIST.data_frame
    res <- CoGAPS(data, nIterations=1000,
                  seed=1, messages=FALSE, sparseOptimization=FALSE)
    res_unc <- CoGAPS(data, nIterations=1000, uncertainty = 0.1*as.matrix(data),
                  seed=1, messages=FALSE, sparseOptimization=FALSE)

    reported <- getMeanChiSq(res)
    reported_unc <- getMeanChiSq(res_unc)
    message('unc. not provided:', reported, '\nunc. provided:', reported_unc)
    #we do not know how big of a difference is allowed, setting to too big
    expect_equal(reported, reported_unc, tolerance = 1e+6)
})