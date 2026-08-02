test_that("featureLoadings and sampleFactors are not all 0s in single-cell", {
  params <- CogapsParams(seed=42,
                         nIterations = 100,
                         nPatterns = 2,
                         sparseOptimization = as.logical(0),
                         distributed="single-cell")
  
  params <- setDistributedParams(params, nSets = 2)
  data(GIST)
  cg <- CoGAPS(GIST.matrix, params=params)

  featureLoadings <- cg@featureLoadings
  sampleFactors <- cg@sampleFactors

  # Check featureLoadings and sampleFactors: smaller dimension is nPatterns
  # and larger dimension matches data dimensions
  expect_false(all(featureLoadings == 0))
  expect_false(all(sampleFactors == 0))
  expect_true(sort((dim(sampleFactors)))[1] == params@nPatterns)
  expect_true(sort((dim(featureLoadings)))[1] == params@nPatterns)
  expect_true(sort((dim(sampleFactors)))[2] == ncol(GIST.matrix))
  expect_true(sort((dim(featureLoadings)))[2] == nrow(GIST.matrix))
})

test_that("featureLoadings and sampleFactors are not all 0s in genome-wide", {
  params <- CogapsParams(seed=42,
                         nIterations = 100,
                         nPatterns = 2,
                         sparseOptimization = as.logical(0),
                         distributed="genome-wide")
  
  params <- setDistributedParams(params, nSets = 2)
  data(GIST)
  cg <- CoGAPS(GIST.matrix, params=params)
  
  featureLoadings <- cg@featureLoadings
  sampleFactors <- cg@sampleFactors

  # Check featureLoadings and sampleFactors: smaller dimension is nPatterns
  # and larger dimension matches data dimensions
  expect_false(all(featureLoadings == 0))
  expect_false(all(sampleFactors == 0))
  expect_true(sort((dim(sampleFactors)))[1] == params@nPatterns)
  expect_true(sort((dim(featureLoadings)))[1] == params@nPatterns)
  expect_true(sort((dim(sampleFactors)))[2] == ncol(GIST.matrix))
  expect_true(sort((dim(featureLoadings)))[2] == nrow(GIST.matrix))
})