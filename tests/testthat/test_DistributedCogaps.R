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

# One consensus pattern makes the stitched matrices one column wide. Without
# drop=FALSE the re-ordering in stitchTogether turned them into vectors, whose
# rownames are NULL, and the run died with "no gene names given".
test_that("a one-pattern distributed run keeps its dim names", {
  data(GIST)
  for (mode in c("genome-wide", "single-cell"))
  {
    params <- CogapsParams(seed=42,
                           nIterations = 100,
                           nPatterns = 1,
                           sparseOptimization = as.logical(0),
                           distributed = mode)
    params <- setDistributedParams(params, nSets = 2)
    cg <- CoGAPS(GIST.matrix, params=params)

    expect_true(is.matrix(cg@featureLoadings))
    expect_true(is.matrix(cg@sampleFactors))
    expect_equal(ncol(cg@featureLoadings), 1L)
    expect_equal(ncol(cg@sampleFactors), 1L)
    expect_equal(nrow(cg@featureLoadings), nrow(GIST.matrix))
    expect_equal(nrow(cg@sampleFactors), ncol(GIST.matrix))
    expect_false(is.null(rownames(cg@featureLoadings)))
    expect_false(is.null(rownames(cg@sampleFactors)))
  }
})

test_that("distributed workers receive subsetted in-memory data", {
  set.seed(1)
  mat <- matrix(rexp(20 * 10), nrow = 20, ncol = 10)

  params <- CogapsParams(seed = 42,
                         nIterations = 30,
                         nPatterns = 2,
                         sparseOptimization = as.logical(0),
                         distributed = "genome-wide")
  params <- setDistributedParams(params, nSets = 2)
  params <- setParam(params, "explicitSets", list(1:10, 11:20))

  observed <- list()
  ns <- asNamespace("CoGAPS")
  original <- get("callInternalCoGAPS", envir = ns)

  wrapper <- function(data, allParams, uncertainty, subsetIndices, workerID,
                     skipInternalSubset=FALSE)
  {
    observed[[length(observed) + 1]] <<- list(
      nrow = nrow(data),
      ncol = ncol(data),
      subsetLength = length(subsetIndices),
      workerID = workerID,
      skipInternalSubset = skipInternalSubset,
      payloadBytes = as.numeric(object.size(data))
    )
    original(data, allParams, uncertainty, subsetIndices, workerID,
             skipInternalSubset)
  }

  assignInNamespace("callInternalCoGAPS", wrapper, ns = "CoGAPS")
  on.exit(assignInNamespace("callInternalCoGAPS", original, ns = "CoGAPS"),
          add = TRUE)

  CoGAPS(mat,
         params = params,
         BPPARAM = BiocParallel::SerialParam(),
         messages = FALSE)

  expect_true(length(observed) > 0)
  expect_true(all(vapply(observed, function(x) x$nrow < nrow(mat), logical(1))))
  expect_true(all(vapply(observed, function(x) x$ncol == ncol(mat), logical(1))))
  expect_true(all(vapply(observed, function(x) x$subsetLength == x$nrow, logical(1))))
  expect_true(all(vapply(observed, function(x) isTRUE(x$skipInternalSubset), logical(1))))

  # benchmark-style check: worker payload should be smaller than full matrix
  fullBytes <- as.numeric(object.size(mat))
  expect_true(all(vapply(observed, function(x) x$payloadBytes < fullBytes, logical(1))))
})

test_that("single-cell distributed workers receive subsetted in-memory data", {
  set.seed(1)
  mat <- matrix(rexp(20 * 10), nrow = 20, ncol = 10)

  params <- CogapsParams(seed = 42,
                         nIterations = 30,
                         nPatterns = 2,
                         sparseOptimization = as.logical(0),
                         distributed = "single-cell")
  params <- setDistributedParams(params, nSets = 2)
  params <- setParam(params, "explicitSets", list(1:5, 6:10))

  observed <- list()
  ns <- asNamespace("CoGAPS")
  original <- get("callInternalCoGAPS", envir = ns)

  wrapper <- function(data, allParams, uncertainty, subsetIndices, workerID,
                     skipInternalSubset=FALSE)
  {
    observed[[length(observed) + 1]] <<- list(
      nrow = nrow(data),
      ncol = ncol(data),
      subsetLength = length(subsetIndices),
      workerID = workerID,
      skipInternalSubset = skipInternalSubset,
      payloadBytes = as.numeric(object.size(data))
    )
    original(data, allParams, uncertainty, subsetIndices, workerID,
             skipInternalSubset)
  }

  assignInNamespace("callInternalCoGAPS", wrapper, ns = "CoGAPS")
  on.exit(assignInNamespace("callInternalCoGAPS", original, ns = "CoGAPS"),
          add = TRUE)

  CoGAPS(mat,
         params = params,
         BPPARAM = BiocParallel::SerialParam(),
         messages = FALSE)

  expect_true(length(observed) > 0)
  expect_true(all(vapply(observed, function(x) x$nrow == nrow(mat), logical(1))))
  expect_true(all(vapply(observed, function(x) x$ncol < ncol(mat), logical(1))))
  expect_true(all(vapply(observed, function(x) x$subsetLength == x$ncol, logical(1))))
  expect_true(all(vapply(observed, function(x) isTRUE(x$skipInternalSubset), logical(1))))

  # benchmark-style check: worker payload should be smaller than full matrix
  fullBytes <- as.numeric(object.size(mat))
  expect_true(all(vapply(observed, function(x) x$payloadBytes < fullBytes, logical(1))))
})
