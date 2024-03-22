test_that("patternMarkers work with threshold = 'all' for general mode", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100,
                  seed=1, messages=FALSE)
    expect_no_error(patternMarkers(res, threshold = "all"))
})

test_that("patternMarkers work with threshold = 'all' for sparse mode", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100,
                  seed=1, messages=FALSE, sparseOptimization=TRUE)
    expect_no_error(patternMarkers(res, threshold = "all"))
    
})

test_that("patternMarkers work with threshold = 'all' and axis = 2", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100,
                  seed=1, messages=FALSE)
    expect_no_error(patternMarkers(res, threshold = "all", axis = 2))
})

test_that("patternMarkers work with threshold = 'all' and axis = 1", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100,
                  seed=1, messages=FALSE)
    expect_no_error(patternMarkers(res, threshold = "all", axis = 1))
})
