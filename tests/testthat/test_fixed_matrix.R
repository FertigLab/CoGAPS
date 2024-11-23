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
  #expect all 0s in the fixed P matrix
  expect_equal(unique(min(res3@sampleFactors), max(res3@sampleFactors)), 0)
})

testthat::test_that("Fixing A matrix works", {
  data("GIST")
  testMatrix <- GIST.matrix
  nPat <- 5
  nIter <- 100
  #run CoGAPS normally
  res1 <- CoGAPS(testMatrix, nPatterns=nPat, nIterations=nIter, seed=42,
                 messages = FALSE)

  #run with fixed A matrix to reconstruct P
  fixedA <- res1@featureLoadings
  param <- CogapsParams()
  param <- setFixedPatterns(param, fixedA, "A")
  res2 <- CoGAPS(testMatrix, param, nPatterns=nPat, nIterations=nIter, seed=42,
                 messages = FALSE)

  #run with random fixedA
  fixedA <- matrix(runif(nrow(testMatrix) * nPat, 1, 10), ncol=nPat)
  param <- CogapsParams()
  param <- setFixedPatterns(param, fixedA, "A")
  res3 <- CoGAPS(testMatrix, param, nPatterns=nPat, nIterations=nIter, seed=42,
                 messages = FALSE)

  #expect that the P with real fixed A is better than P with random fixed A
  expect_true(abs(sum((res1@sampleFactors - res2@sampleFactors))) <
              abs(sum((res1@sampleFactors - res3@sampleFactors)))
  )
  #expect all 0s in the fixed A matrix
  expect_equal(unique(min(res3@featureLoadings), max(res3@featureLoadings)), 0)
})

test_that("chisq changes between iterations with fixed A", {
  data("GIST")
  testMatrix <- GIST.matrix
  nPat <- 5
  nIter <- 100
  #run CoGAPS normally
  res1 <- CoGAPS(testMatrix, nPatterns=nPat, nIterations=nIter, seed=42,
                 outputFrequency = 10, messages = FALSE)

  #run with fixed A
  fixed <- res1@featureLoadings
  param <- CogapsParams()
  param <- setFixedPatterns(param, fixed, "A")
  res2 <- CoGAPS(testMatrix, param, nPatterns=nPat, nIterations=nIter, seed=42,
                 outputFrequency = 10, messages = FALSE)

  expect_true(length(unique(res2@metadata$chisq))==length(res2@metadata$chisq))
})

test_that("chisq changes between iterations with fixed P", {
  data("GIST")
  testMatrix <- GIST.matrix
  nPat <- 5
  nIter <- 100
  #run CoGAPS normally
  res1 <- CoGAPS(testMatrix, nPatterns=nPat, nIterations=nIter, seed=42,
                 outputFrequency = 10, messages = FALSE)

  #run with fixed P
  fixed <- res1@sampleFactors
  param <- CogapsParams()
  param <- setFixedPatterns(param, fixed, "P")
  res2 <- CoGAPS(testMatrix, param, nPatterns=nPat, nIterations=nIter, seed=42,
                 outputFrequency = 10, messages = FALSE)

  expect_true(length(unique(res2@metadata$chisq))==length(res2@metadata$chisq))
})
