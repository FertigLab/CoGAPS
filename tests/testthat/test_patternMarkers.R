
test_that("patternMarkers work with threshold = 'all' for sparse mode", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100,
                  seed=1, messages=FALSE, sparseOptimization=TRUE)
    expect_no_error(patternMarkers(res, threshold = "all"))
    
})

test_that("patternMarkers work with threshold = 'all'", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100,
                  seed=1, messages=FALSE)
    expect_no_error(patternMarkers(res, threshold = "all"))
})

test_that("test outputs generic in threshold = all", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100,
                  seed=1, messages=FALSE)
    test <- patternMarkers(res, threshold = "all")
    marker_lengths <- sapply(test$PatternMarkers, length)
    #not all patterns have no markers
    expect_equal(sum(marker_lengths > 0), length(marker_lengths))
    #all markers are non-empty
    expect_equal(sum(nchar(names(test$PatternMarkers)) > 0),
                 length(names(test$PatternMarkers)))
    expect_equal(length(unlist(test$PatternMarkers)),
                 nrow(res@featureLoadings))
})


#functional tests
gapsMock <- function(mock){
    obj <- new(
    "CogapsResult",
    Amean = mock$featureLoadings,
    Pmean = mock$sampleFactors,
    Asd = mock$featureLoadings, # just putting in, these values arent used
    Psd = mock$sampleFactors, # just putting in, these values arent used
    meanChiSq = mock$meanChiSq,
    geneNames = mock$geneNames,
    sampleNames = mock$sampleNames,
    diagnostics = NULL
  )
  obj
}

test_that("patternMarkers work with threshold = 'all' for mock object", {
  # mock CogapsResult object for functional tests
  mock <- list(
    featureLoadings = diag(1, 5, 5), #make each nth gene marker in nth pattern
    sampleFactors = matrix(rep(1, 25), nrow = 5), # 5x5 matrix of 1s
    factorStdDev = matrix(runif(25), nrow = 5), # 5x5 matrix of random numbers
    meanChiSq = runif(1), # single random number
    geneNames = paste0("Gene", 1:5), # vector of gene names
    sampleNames = paste0("Sample", 1:5) # vector of sample names
  )

  # create a new CogapsResult object
  obj <- gapsMock(mock)

  pm <- patternMarkers(obj, threshold = "all")
  expect_equal(pm$PatternMarkers$Pattern_1, "Gene1")
  expect_equal(pm$PatternMarkers$Pattern_2, "Gene2")
  expect_equal(pm$PatternMarkers$Pattern_3, "Gene3")
  expect_equal(pm$PatternMarkers$Pattern_4, "Gene4")
  expect_equal(pm$PatternMarkers$Pattern_5, "Gene5")
})

test_that("genes are sorted by by their ranks in the output", {
  ## a more complicated mock object
  test_mat <- t(matrix(c(c(10,1,1,1,1),
                        c(2,10,2,2,2),
                        c(3,3,10,3,3),
                        c(4,4,4,10,4),
                        c(5,5,5,5,10)), nrow=5, ncol=5))
  test_mat <- rbind(test_mat,diag(1, 5, 5)[,5:1])
  mock <- list(
    featureLoadings = test_mat,
    sampleFactors = matrix(rep(1, 25), nrow = 5), # 5x5 matrix of 1s
    factorStdDev = matrix(runif(25), nrow = 5), # 5x5 matrix of random numbers
    meanChiSq = runif(1), # single random number
    geneNames = paste0("Gene", 1:10), # vector of gene names
    sampleNames = paste0("Sample", 1:5) # vector of sample names
  )
  obj <- gapsMock(mock)

  pm <- patternMarkers(obj, threshold = "all")

  expect_true(all.equal(pm$PatternMarkers$Pattern_1, c("Gene10","Gene1")))
  expect_true(all.equal(pm$PatternMarkers$Pattern_3, c("Gene8","Gene3")))
  expect_true(all.equal(pm$PatternMarkers$Pattern_5, c("Gene6","Gene5")))

})

test_that("all patterns are returned even if no markers",{
  mock <- list(
    featureLoadings = diag(1, 5, 5)[,5:1], #make each nth gene marker in nth pattern
    sampleFactors = matrix(rep(1, 25), nrow = 5), # 5x5 matrix of 1s
    factorStdDev = matrix(runif(25), nrow = 5), # 5x5 matrix of random numbers
    meanChiSq = runif(1), # single random number
    geneNames = paste0("Gene", 1:5), # vector of gene names
    sampleNames = paste0("Sample", 1:5) # vector of sample names
  )
  #make gene1 not marker of any pattern
  mock$featureLoadings[1:2,] <- c(0,0,0,0,0)

  obj <- gapsMock(mock)
  pm <- patternMarkers(obj, threshold = "all")

  expect_true(length(pm$PatternMarkers) == 5)
})

test_that("patternMarkers works with lp", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100, nPatterns=5,
                  seed=1, messages=FALSE, sparseOptimization=TRUE)
    expect_no_error(patternMarkers(res, lp=list(my_lp=c(1,0,0,0,0))))
    
})