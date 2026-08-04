context("CoGAPS")

# Accessors and analysis methods on CogapsResult that the suite did not touch:
# getFeatureLoadings, getVersion, getOriginalParameters, calcZ, reconstructGene,
# binaryA, calcCoGAPSStat, toCSV/fromCSV.

data(GIST)
# 1000 iterations: with fewer, part of the standard-deviation matrix is still
# exactly zero and calcZ warns about it, which would make these tests noisy.
res <- CoGAPS(GIST.matrix, nPatterns=3, nIterations=1000, seed=1, messages=FALSE)

test_that("the matrix accessors return the corresponding slots",
{
    expect_identical(getFeatureLoadings(res), res@featureLoadings)
    expect_identical(getSampleFactors(res), res@sampleFactors)
    expect_identical(getAmplitudeMatrix(res), res@featureLoadings)
    expect_identical(getPatternMatrix(res), res@sampleFactors)

    expect_equal(nrow(getFeatureLoadings(res)), nrow(GIST.matrix))
    expect_equal(ncol(getFeatureLoadings(res)), 3)
})

test_that("getVersion and getOriginalParameters describe the run",
{
    # getVersion returns the package_version recorded in the result metadata
    expect_true(is(getVersion(res), "numeric_version"))
    expect_equal(format(getVersion(res)), as.character(packageVersion("CoGAPS")))
    params <- getOriginalParameters(res)
    expect_true(is(params, "CogapsParams"))
    expect_equal(params@nPatterns, 3)
    expect_equal(params@seed, 1)
})

test_that("calcZ returns mean/stddev z-scores of the requested matrix",
{
    zA <- calcZ(res, "featureLoadings")
    zP <- calcZ(res, "sampleFactors")
    expect_equal(dim(zA), dim(res@featureLoadings))
    expect_equal(dim(zP), dim(res@sampleFactors))
    expect_true(all(is.finite(zA)))
    expect_true(all(is.finite(zP)))
    expect_error(calcZ(res, "notAMatrix"), "whichMatrix")
})

test_that("calcZ warns when the standard deviation matrix contains zeros",
{
    # short runs leave exact zeros in the sd matrix, which would divide by zero
    zeroed <- res
    zeroed@loadingStdDev[1, 1] <- 0
    expect_warning(calcZ(zeroed, "featureLoadings"), "standard deviation")
})

test_that("reconstructGene rebuilds rows of the data matrix",
{
    full <- reconstructGene(res)
    expect_equal(dim(full), dim(GIST.matrix))
    # the reconstruction is A %*% t(P)
    expect_equal(unname(as.matrix(full)),
                 unname(res@featureLoadings %*% t(res@sampleFactors)),
                 tolerance=1e-4)

    one <- reconstructGene(res, genes=1:5)
    expect_equal(nrow(one), 5)
})

test_that("binaryA draws the thresholded amplitude heatmap",
{
    # binaryA is a plotting function -- it thresholds calcZ() of the A matrix and
    # draws a heatmap, returning whatever mtext() returns rather than the matrix.
    # Regression: it called calcZ(object) without the mandatory whichMatrix, so
    # any call failed with 'argument "whichMatrix" is missing'.
    # draw into a throwaway file so the run leaves no Rplots.pdf behind
    tmp <- tempfile(fileext=".pdf")
    pdf(tmp)
    on.exit({ dev.off(); unlink(tmp) }, add=TRUE)
    expect_error(binaryA(res, threshold=3), NA)
})

test_that("calcCoGAPSStat scores a gene set",
{
    genes <- rownames(GIST.matrix)
    sets <- list(setA=genes[1:50], setB=genes[51:100])
    stat <- calcCoGAPSStat(res, sets=sets, whichMatrix="featureLoadings",
                           numPerm=100)
    expect_true(is.list(stat) || is.data.frame(stat))
    expect_true(length(stat) > 0)
})

test_that("toCSV and fromCSV round-trip a result",
{
    dir <- file.path(tempdir(), "cogaps_csv_roundtrip")
    dir.create(dir, showWarnings=FALSE)
    on.exit(unlink(dir, recursive=TRUE), add=TRUE)

    toCSV(res, dir)
    expect_true(all(file.exists(file.path(dir,
        c("featureLoadings.csv", "sampleFactors.csv",
          "loadingStdDev.csv", "factorStdDev.csv")))))

    back <- fromCSV(dir)
    expect_true(is(back, "CogapsResult"))
    expect_equal(dim(back@featureLoadings), dim(res@featureLoadings))
    expect_equal(dim(back@sampleFactors), dim(res@sampleFactors))

    # the round-trip must preserve the matrix type of the slots, not hand back
    # the data.frame read.csv produces
    expect_true(is.matrix(back@featureLoadings))
    expect_true(is.matrix(back@sampleFactors))
    expect_equal(unname(back@featureLoadings), unname(res@featureLoadings),
                 tolerance=1e-6)
    expect_equal(unname(back@sampleFactors), unname(res@sampleFactors),
                 tolerance=1e-6)
})

test_that("show() on CogapsParams works with checkpoint parameters set",
{
    # regression: the checkpointInFile branch of show() referenced a bare
    # `checkpointInFile` instead of `object@checkpointInFile`, so printing any
    # params object with a checkpoint file set raised "object not found"
    p <- CogapsParams(nPatterns=3)
    p <- setParam(p, "checkpointInFile", "somewhere.out")
    expect_output(show(p), "checkpointInFile")

    p2 <- CogapsParams(nPatterns=3)
    p2 <- setParam(p2, "checkpointOutFile", "out.out")
    expect_output(show(p2), "checkpointOutFile")
})
