test_that("getPatternGeneSet works on enrichment test", {
  #set up
  data(GIST)
  gs.test <- list(
    "gs1" = c("Hs.2", "Hs.4", "Hs.36", "Hs.96", "Hs.202"),
    "gs2" = c("Hs.699463", "Hs.699288", "Hs.699280", "Hs.699154", "Hs.697294")
  )
  expect_no_error(getPatternGeneSet(object = GIST.result, gene.sets = gs.test, method = "enrichment"))
})

test_that("getPatternGeneSet works on overrepresentation test", {
  #set up
  data(GIST)
  gs.test <- list(
    "gs1" = c("Hs.2", "Hs.4", "Hs.36", "Hs.96", "Hs.202"),
    "gs2" = c("Hs.699463", "Hs.699288", "Hs.699280", "Hs.699154", "Hs.697294")
  )
  expect_no_error(getPatternGeneSet(object = GIST.result, gene.sets = gs.test, method = "overrepresentation"))
})

test_that("plotPatternGeneSet renders a plot for enrichment test", {
  #set up
  data(GIST)
  gs.test <- list(
    "gs1" = c("Hs.2", "Hs.4", "Hs.36", "Hs.96", "Hs.202"),
    "gs2" = c("Hs.699463", "Hs.699288", "Hs.699280", "Hs.699154", "Hs.697294")
  )
  gpgs_res <- getPatternGeneSet(object = GIST.result, gene.sets = gs.test, method = "enrichment")
  pl <- plotPatternGeneSet(patterngeneset = gpgs_res, whichpattern = 1, padj_threshold = 1)
  
  expect_is(pl$layers[[1]], "ggproto")
  expect_match(pl$labels$title, "Enriched gene sets in Pattern_1")
})

test_that("plotPatternGeneSet renders a plot for overrepresentation test", {
  #set up
  data(GIST)
  gs.test <- list(
    "gs1" = c("Hs.2", "Hs.4", "Hs.36", "Hs.96", "Hs.202"),
    "gs2" = c("Hs.699463", "Hs.699288", "Hs.699280", "Hs.699154", "Hs.697294")
  )
  gpgs_res <- getPatternGeneSet(object = GIST.result, gene.sets = gs.test, method = "overrepresentation", threshold = "cut")
  pl <- plotPatternGeneSet(patterngeneset = gpgs_res, whichpattern = 1, padj_threshold = 1)
  
  expect_is(pl$layers[[1]], "ggproto")
  expect_match(pl$labels$title, "Overrepresented gene sets in Pattern_1")
})
