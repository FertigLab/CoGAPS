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
a custom (as compared to `testthat::use_catch()`) runner to run  cpp tests, exposes `run_catch_unit_tests()` in R, usage: 
```
#just run all cpp tests
CoGAPS::run_catch_unit_tests()

#use cpp xml reporter to see debug info
CoGAPS:::run_catch_unit_tests(reporter=“xml")

#call a single cpp test by name
CoGAPS:::run_catch_unit_tests_by_tag("Test Vector.h")

#call cpp test(s) by tag "vector"
CoGAPS:::run_catch_unit_tests_by_tag("[vector]")

```
The tags need to be defined in `/src/cpp_tests/[test-name].cpp`

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

After changes are done, update `configure.ac` by running 

in terminal:
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


