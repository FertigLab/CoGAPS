
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

test_that("no empty patternMarkers and their names in threshold = all", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100,
                  seed=1, messages=FALSE)
    test <- patternMarkers(res, threshold = "all")
    marker_lengths <- sapply(test$PatternMarkers, length)
    expect_equal(sum(marker_lengths > 0), length(marker_lengths))
    expect_equal(sum(nchar(names(test$PatternMarkers)) > 0),
                 length(names(test$PatternMarkers)))
})

test_that("no empty patternMarkers and their names in threshold = cut", {
    #set up
    data(GIST)
    res <- CoGAPS(GIST.data_frame, nIterations=100,
                  seed=1, messages=FALSE)
    test <- patternMarkers(res, threshold = "cut")
    marker_lengths <- sapply(test$PatternMarkers, length)
    expect_equal(sum(marker_lengths > 0), length(marker_lengths))
    expect_equal(length(test$PatternMarkers),
                 length(names(test$PatternMarkers)))
})

# mock CogapsResult object for functional tests
# a list that will become CogapsResult object
mock <- list(
  featureLoadings = diag(1, 5, 5), #make each nth gene marker in nth pattern
  sampleFactors = matrix(rep(1, 25), nrow = 5), # 5x5 matrix of 1s
  factorStdDev = matrix(runif(25), nrow = 5), # 5x5 matrix of random numbers
  meanChiSq = runif(1), # single random number
  geneNames = paste0("Gene", 1:5), # vector of gene names
  sampleNames = paste0("Sample", 1:5) # vector of sample names
)

# create a new CogapsResult object
obj <- new(
  "CogapsResult",
  Amean = mock$featureLoadings,
  Pmean = mock$sampleFactors,
  Asd = mock$factorStdDev,
  Psd = mock$factorStdDev,
  meanChiSq = mock$meanChiSq,
  geneNames = mock$geneNames,
  sampleNames = mock$sampleNames,
  diagnostics = NULL
)

test_that("patternMarkers work with threshold = 'cut' for mock object", {
  pm <- patternMarkers(obj, threshold = "cut")
  expect_equal(pm$PatternMarkers$Pattern_1, "Gene1")
  expect_equal(pm$PatternMarkers$Pattern_2, "Gene2")
  expect_equal(pm$PatternMarkers$Pattern_3, "Gene3")
  expect_equal(pm$PatternMarkers$Pattern_4, "Gene4")
  expect_equal(pm$PatternMarkers$Pattern_5, "Gene5")
})


test_that("patternMarkers work with threshold = 'all' for mock object", {
  pm <- patternMarkers(obj, threshold = "all")
  expect_equal(pm$PatternMarkers$Pattern_1, "Gene1")
  expect_equal(pm$PatternMarkers$Pattern_2, "Gene2")
  expect_equal(pm$PatternMarkers$Pattern_3, "Gene3")
  expect_equal(pm$PatternMarkers$Pattern_4, "Gene4")
  expect_equal(pm$PatternMarkers$Pattern_5, "Gene5")
})
