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
  
  # Check that featureLoadings and sampleFactors are not all zeros
  expect_false(all(featureLoadings == 0))
  expect_false(all(sampleFactors == 0))
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
  
  # Check that featureLoadings and sampleFactors are not all zeros
  expect_false(all(featureLoadings == 0))
  expect_false(all(sampleFactors == 0))
})