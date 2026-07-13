context("CoGAPS")

# Recompute the mean chi-square directly from the returned mean factor matrices
# and the default uncertainty model, and check it against the value CoGAPS
# reports (GapsStatistics::meanChiSq).
#
# The reconstruction is  M = A %*% t(P)  where A = featureLoadings (mean of the A
# matrix) and P = sampleFactors (mean of the P matrix). The uncertainty is the
# relative-error model floored at 0.1:  S = max(0.1 * D, 0.1). This model is the
# same for the dense and sparse samplers (see DenseNormalModel / SparseNormalModel
# and issues #17/#19), so the identical recomputation must hold in both modes.
manualMeanChiSq <- function(res, D)
{
    A <- res@featureLoadings          # genes x patterns   (mean of A)
    P <- res@sampleFactors            # samples x patterns  (mean of P)
    S <- pmax(0.1 * D, 0.1)
    sum(((D - A %*% t(P)) / S)^2)
}

test_that("chi-square reported by CoGAPS matches manually calculated (w/uncertainty)",
{
    data(GIST)
    D <- as.matrix(GIST.matrix)

    # dense sampler
    res <- CoGAPS(GIST.matrix, nPatterns=5, nIterations=500, outputFrequency=100,
        seed=42, messages=FALSE)
    expect_equal(getMeanChiSq(res), manualMeanChiSq(res, D), tolerance=1e-3)

    # sparse sampler -- same uncertainty model, so the same recomputation applies
    res_sp <- CoGAPS(GIST.matrix, nPatterns=5, nIterations=500, outputFrequency=100,
        seed=42, messages=FALSE, sparseOptimization=TRUE)
    expect_equal(getMeanChiSq(res_sp), manualMeanChiSq(res_sp, D), tolerance=1e-3)
})
