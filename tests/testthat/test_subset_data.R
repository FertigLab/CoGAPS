context("CoGAPS")

test_that("standard cogaps on a subset of the data",
{
    data(GIST)    
    subset <- sample(1:nrow(GIST.matrix), 500)
    result <- CoGAPS(GIST.matrix, nIterations=50, messages=FALSE, seed=42,
        subsetIndices=subset, subsetDim=1)
    expect_equal(length(subset), nrow(result@featureLoadings))
})

test_that("subsetting data with explicit sets",
{
    gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")

    # distributed cogaps across features
    in_sets <- list(1:225, 226:450, 451:675, 676:900)
    result <- CoGAPS(gistMtxPath, nPatterns=3, explicitSets=in_sets,
        nIterations=200, messages=FALSE, seed=42, distributed="genome-wide")
    featureNames <- rownames(result@featureLoadings)
    out_sets <- lapply(getSubsets(result), function(set) which(featureNames %in% set))
    expect_true(all(sapply(1:4, function(i) all.equal(out_sets[[i]], in_sets[[i]]))))

    # distributed cogaps across samples
    in_sets <- list(1:225, 226:450, 451:675, 676:900)
    result <- CoGAPS(gistMtxPath, nPatterns=3, explicitSets=in_sets, seed=42,
        nIterations=200, messages=FALSE, distributed="single-cell",
        transposeData=TRUE, singleCell=TRUE)
    sampleNames <- rownames(result@sampleFactors)
    out_sets <- lapply(getSubsets(result), function(set) which(sampleNames %in% set))
    expect_true(all(sapply(1:4, function(i) all.equal(out_sets[[i]], in_sets[[i]]))))
})

test_that("subsetting data with uniform sets",
{
    gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")

    # distributed cogaps across features
    result <- CoGAPS(gistMtxPath, nPatterns=3, nIterations=200, messages=FALSE,
        seed=42, distributed="genome-wide")
    featureNames <- rownames(result@featureLoadings)
    sets <- lapply(getSubsets(result), function(set) which(featureNames %in% set))
    expect_equal(sum(sapply(sets, length)), nrow(result@featureLoadings))

    # distributed cogaps across samples
    result <- CoGAPS(gistMtxPath, nPatterns=3, nIterations=200, messages=FALSE,
        seed=42, distributed="single-cell", transposeData=TRUE, singleCell=TRUE)
    sampleNames <- rownames(result@sampleFactors)
    sets <- lapply(getSubsets(result), function(set) which(sampleNames %in% set))
    expect_equal(sum(sapply(sets, length)), nrow(result@sampleFactors))
})

test_that("subsetting data with annotation weights",
{
    # TODO address how weighted sampling works with duplicates, do we need to
    # allow passing a value for setSize in this case?
    # we should collapse down using the mean
    # prevent multiple copies from being in the same set

    #data(GIST)
    #gistMtxPath <- system.file("extdata/GIST.mtx", package="CoGAPS")    
#
    ## create annotations
    #weight <- c(1, 2, 3)
    #names(weight) <- c("A", "B", "C")
    #anno <- sample(names(weight), nrow(GIST.matrix), replace=TRUE)
    #params <- CogapsParams()
    #params <- setAnnotationWeights(params, anno, weight)
#
    ## distributed cogaps across features
    #result <- CoGAPS(gistMtxPath, params, nPatterns=3, nIterations=200,
    #    messages=TRUE, seed=42, distributed="genome-wide")
    #featureNames <- rownames(result@featureLoadings)
    #sets <- lapply(getSubsets(result), function(set) which(featureNames %in% set))
    #expect_equal(nrow(result@featureLoadings), nrow(GIST.matrix))
    #expect_equal(sum(sapply(sets, length)), nrow(result@featureLoadings))
#
    ## distributed cogaps across samples
    #result <- CoGAPS(gistMtxPath, params, nPatterns=3, nIterations=200,
    #    messages=FALSE, seed=42, distributed="single-cell", transposeData=TRUE,
    #    singleCell=TRUE)
    #sampleNames <- rownames(result@sampleFactors)
    #sets <- lapply(getSubsets(result), function(set) which(sampleNames %in% set))
    #expect_equal(nrow(result@sampleFactors), nrow(GIST.matrix))
    #expect_equal(sum(sapply(sets, length)), nrow(result@sampleFactors))
})