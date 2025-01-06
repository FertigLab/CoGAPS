context("CoGAPS")

test_that("CogapsParams class",
{
  params<-new("CogapsParams")
  cat("\n")
  print(params)
  expect_true("CogapsParams" %in% class(params))
})