context("CoGAPS")

test_that("Npatterns are required input",
{
    data(GIST)
    expect_error(CoGAPS(data=GIST.data_frame, nIterations=100, messages=FALSE), "nPatterns")
    expect_no_error(CoGAPS(data=GIST.data_frame, nPatterns=7, nIterations=100, messages=FALSE))
})

test_that("CogapsParams class",
{
  params<-new("CogapsParams", nPatterns=7)
  cat("\n")
  print(params)
  expect_true("CogapsParams" %in% class(params))
})