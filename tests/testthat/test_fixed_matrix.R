test_that("Fixing P matrix works", {
  data("GIST")
  testMatrix <- GIST.matrix
  nPat <- 5
  nIter <- 100
  #run CoGAPS normally
  res1 <- CoGAPS(testMatrix, nPatterns=nPat, nIterations=nIter, seed=42,
                 messages = FALSE)

  #run with fixed P matrix to reconstruct A
  fixedP <- res1@sampleFactors
  param <- CogapsParams()
  param <- setFixedPatterns(param, fixedP, "P")
  res2 <- CoGAPS(testMatrix, param, nPatterns=nPat, nIterations=nIter, seed=42,
                 messages = FALSE)

  #run with random fixedP
  fixedP <- matrix(runif(ncol(testMatrix) * nPat, 1, 10), ncol=nPat)
  param <- CogapsParams()
  param <- setFixedPatterns(param, fixedP, "P")
  res3 <- CoGAPS(testMatrix, param, nPatterns=nPat, nIterations=nIter, seed=42,
                 messages = FALSE)

  #expect that the A with real fixed P is better than A with random fixed P
  expect_true(abs(sum((res1@featureLoadings - res2@featureLoadings))) <
              abs(sum((res1@featureLoadings - res3@featureLoadings))))
})