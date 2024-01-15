test_that("patternMarkers work with threshold = 'all'", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100, outputFrequency=50, seed=1, messages=FALSE)

    test <- patternMarkers(res, threshold = "all")
    expect_gt(length(test$PatternMarkers), 0)
})