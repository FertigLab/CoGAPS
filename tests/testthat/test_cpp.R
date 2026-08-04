context("C++")

# The Catch2 C++ suite is run from here so that a C++ failure fails
# `R CMD check`. The report is written as XML to a temporary file and parsed, so
# that every C++ TEST_CASE shows up as its own testthat expectation instead of
# the whole suite collapsing into a single pass/fail. Writing to a file also
# keeps Catch's output out of stdout -- it is emitted by C++, so testthat cannot
# capture it.
#
# To run the same tests interactively, with the ordinary human-readable output,
# use the console form (see src/README.md):
#
#     CoGAPS:::run_catch_unit_tests()
#     CoGAPS:::run_catch_unit_tests_by_tag("[vector]")

test_that("the C++ unit tests are compiled into this build",
{
    cases <- CoGAPS:::catch_test_case_names()
    if (length(cases) == 0)
        skip(paste("no C++ test cases are registered -- this build was made",
                   "with --disable-cpp-tests, or on Windows, where",
                   "src/Makevars.win lists no cpp_tests objects"))
    expect_gt(length(cases), 0)
})

test_that("Catch unit tests pass",
{
    cases <- CoGAPS:::catch_test_case_names()
    if (length(cases) == 0)
        skip("no C++ test cases are registered")

    # the C++ file-parser tests read these data paths from the global environment
    gistCsvPath <<- system.file("extdata/GIST.csv", package="CoGAPS")
    gistTsvPath <<- system.file("extdata/GIST.tsv", package="CoGAPS")
    gistMtxPath <<- system.file("extdata/GIST.mtx", package="CoGAPS")
    gistGctPath <<- system.file("extdata/GIST.gct", package="CoGAPS")

    reportFile <- tempfile(fileext=".xml")
    on.exit(unlink(reportFile), add=TRUE)

    # returns the number of failed assertions (capped at 255); 0 == all pass
    nFailed <- CoGAPS:::run_catch_unit_tests(reporter="xml", output=reportFile)

    # xml2 is only a suggested dependency, so fall back to the summary check
    if (!requireNamespace("xml2", quietly=TRUE))
    {
        expect_equal(nFailed, 0L)
        return(invisible(NULL))
    }

    report <- xml2::read_xml(reportFile)
    testCases <- xml2::xml_find_all(report, "//TestCase")
    expect_equal(length(testCases), length(cases))

    for (tc in testCases)
    {
        name <- xml2::xml_attr(tc, "name")
        # OverallResult (singular) is the per-case verdict; OverallResults
        # (plural) is the assertion tally on sections and on the whole group
        verdict <- xml2::xml_find_first(tc, "./OverallResult")
        passed <- identical(xml2::xml_attr(verdict, "success"), "true")

        if (!passed)
        {
            where <- paste0(xml2::xml_attr(tc, "filename"), ":",
                            xml2::xml_attr(tc, "line"))
            fail(paste0("C++ test case failed: ", name, " (", where, ")"))
        }
        else
        {
            succeed()
        }
    }

    # belt and braces: the tally has to agree with the per-case verdicts
    expect_equal(nFailed, 0L)
})
