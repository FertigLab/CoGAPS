testthat('chi-square reported by CoGAPS mathches manually calculated one',{
    data(GIST)
    data <- GIST.data_frame
    res <- CoGAPS(data, nIterations=1000,
                  seed=1, messages=FALSE, sparseOptimization=FALSE)
    reported <- getMeanChiSq(res)

    A <- getAmplitudeMatrix(res)
    P <- getPatternMatrix(res)
    M <- A %*% t(P)

    calculated <- sum((M - data)^2 / data)

    message("Reported:", reported, "\nCalculated:", calculated)

    expect_equal(reported, calculated, tolerance = 1e-6)
})
