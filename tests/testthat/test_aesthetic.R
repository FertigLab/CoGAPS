context("lints")

test_that("Package Style", {
    lintr::expect_lint_free()
})
