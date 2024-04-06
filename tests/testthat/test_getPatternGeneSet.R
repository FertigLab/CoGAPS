test_that("getPatternGeneSet works on enrichment test", {
  #set up
  data(GIST)
  gs.test <- list(
    "gs1" = c("Hs.2", "Hs.4", "Hs.36", "Hs.96", "Hs.202"),
    "gs2" = c("Hs.699463", "Hs.699288", "Hs.699280", "Hs.699154", "Hs.697294")
  )
  expect_error(getPatternGeneSet(object = GIST.result, gene.sets = gs.test, method = "enrichment"), NA)
})

test_that("getPatternGeneSet works on overrepresentation test", {
  #set up
  data(GIST)
  gs.test <- list(
    "gs1" = c("Hs.2", "Hs.4", "Hs.36", "Hs.96", "Hs.202"),
    "gs2" = c("Hs.699463", "Hs.699288", "Hs.699280", "Hs.699154", "Hs.697294")
  )
  expect_error(getPatternGeneSet(object = GIST.result, gene.sets = gs.test, method = "overrepresentation"), NA)
})
