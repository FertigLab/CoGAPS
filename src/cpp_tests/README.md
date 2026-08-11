# CPP tests in CoGAPS

## cpp unit tests
The cpp unit tests are located in the `/src/cpp_tests`. As at 2024-10-02, the tests are oudated and need refactoring. A separate project will be launched to refactor the cpp tests.

## launching cpp tests from R
In `#7d00008` a way to launch individual cpp unit tests has been introduced by @favorov, below is he list of the most important components that make it possible. The cpp tests use `catch` that ships with `testthat` R package.


## testthat-tweak.h
makes using cpp tags in testthat (context-style) possible,
usage: add the following to each individual cpp test.
```
include testthat.h
include testthat-tweak.h
```
the cpp tests that do not have SECTION will not be called
by testthat


## test-runner.cpp
a custom (as compared to `testthat::use_catch()`) runner to run cpp tests. It
exposes three functions in R. They are internal, so reach them with `:::`.

### Running them by hand -- the plain, readable way

This is what you want while working on the C++ code. Output goes to the console
in the usual Catch format; the return value is the number of failed assertions,
0 meaning everything passed.

```
#just run all cpp tests
CoGAPS:::run_catch_unit_tests()

#use cpp xml reporter to see debug info
CoGAPS:::run_catch_unit_tests(reporter="xml")

#call a single cpp test by name
CoGAPS:::run_catch_unit_tests_by_tag("Test Vector.h")

#call cpp test(s) by tag "vector"
CoGAPS:::run_catch_unit_tests_by_tag("[vector]")

#call cpp test(s) by tag "vector" or "green"
CoGAPS:::run_catch_unit_tests_by_tag("[vector],[green]")

#call cpp test(s) by tag "vector" and "green"
CoGAPS:::run_catch_unit_tests_by_tag("[vector][green]")
CoGAPS:::run_catch_unit_tests_by_tag("[green][vector]")

```

Some of the cpp tests read the paths of the packaged data files from the global
environment, so set those first or the file-parser cases will fail:

```
gistCsvPath <<- system.file("extdata/GIST.csv", package="CoGAPS")
gistTsvPath <<- system.file("extdata/GIST.tsv", package="CoGAPS")
gistMtxPath <<- system.file("extdata/GIST.mtx", package="CoGAPS")
gistGctPath <<- system.file("extdata/GIST.gct", package="CoGAPS")
```

### How testthat runs them

`tests/testthat/test_cpp.R` runs the same suite as part of `devtools::test()` and
`R CMD check`, so a broken C++ test breaks the R build. It does not use the
console form: Catch writes its report from C++, straight to stdout, where
testthat cannot capture it, and the whole suite would collapse into one
pass/fail. Instead it asks for the xml reporter and a file:

```
reportFile <- tempfile(fileext=".xml")
CoGAPS:::run_catch_unit_tests(reporter="xml", output=reportFile)
```

`output=""` (the default) keeps writing to stdout, which is why the plain call
above still behaves the way it always has. The test then parses the file with
`xml2` and turns every `<TestCase>` into its own testthat expectation, so a
failure is reported by name and source location rather than as "expected 0, got
N".

`catch_test_case_names()` returns the names of every TEST_CASE compiled in. The
test uses it to tell "all C++ tests passed" apart from "no C++ tests were built"
-- with `--disable-cpp-tests`, or on Windows where `Makevars.win` lists no
`cpp_tests` objects, the suite is empty and would otherwise pass vacuously. In
that case the test skips instead, with an explanation.

```
length(CoGAPS:::catch_test_case_names())   # 58 as of this writing
```

The tags need to be defined in `/src/cpp_tests/[test-name].cpp` in TEST_CASE:

```
TEST_CASE("Mt test case -- what is it","[tag1][tag2]")
{
}
```

## adding cpp tests to compilation / making changes

To ask for the cpp tests to be compiled, each `test.cpp` needs to be added to the `configure.ac`. Example from `configure.ac`:

```
# add c++ tests to source list
if test "x$cpp_tests_disable" != "xyes" ; then
    echo "Enabling C++ Unit Tests"
    GAPS_CPP_FLAGS+=" -DGAPS_CPP_UNIT_TESTS "
    GAPS_SOURCE_FILES+=" cpp_tests/testVector.o"
fi
```

In the above example, the `testVector.o` is being added and will be compiled (if test compilation is not disabled).

After changes are done, update configure from `configure.ac` by running autoconf in terminal:

```
autoconf
```

in R session:

```
Rcpp::compileAttributes()
```

## disabling cpp test compilation
It may be needed to disable the compilation of tests. There is a specific compilation parameter that controls it. Again, this parameter is set in `configure.ac` (see above). 
Usage:

```
#install with tests disabled
options(configure.args = list(CoGAPS = "--disable-cpp-tests"))

#another way
devtools::install(args = "--configure-args='--enable-cpp-tests'")
```


