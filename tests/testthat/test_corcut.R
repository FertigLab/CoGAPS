test_that("corcut splits patterns according to cut", {
  opt <- expand.grid(n=c(2, 5, 10), len=c(10, 100),
                     s=c(2, 5, 10), jit = c(0.1, 0.5, 1))
  for (i in seq_along(nrow(opt))) {
    sets <- CoGAPS:::makeCorPatternSets(n=opt$n[i], len=opt$len[i],
                                        s=opt$s[i], jit=opt$jit[i])
    sets <- do.call(cbind, sets)
    clusters <- CoGAPS:::corcut(sets, opt$n[i], minNS=1)
    expect_equal(length(clusters), opt$n[i])
  }
})

test_that("corcut correctly handles a dataset with NA values", {
  patterns <- matrix(c(1, 2, NA, 4, 5, 6), nrow=2)
  expect_error(CoGAPS:::corcut(patterns, cut=2, minNS=1),
               "NA values in correlation of patterns")
})

test_that("corcut correctly handles a dataset with negative correlations", {
  patterns <- matrix(c(1, 2, 3, -4, -5, -6), nrow=2)
  clusters <- CoGAPS:::corcut(patterns, cut=2, minNS=1)
  expect_equal(length(clusters), 2)
})