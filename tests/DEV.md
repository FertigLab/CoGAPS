## testthat-tweak.h
makes using cpp tags in testthat (context-style) possible,
used intead of `include catch.h`
usage: 
```
include testthat.h
include testthat-tweak.h
```

the cpp tests that do not have SECTION will not be called
by testthat


## test-runner.cpp
a custom (as compared to use-catch)

to run all cpp tests - usage: 
```
CoGAPS::run_catch_unit_tests()
```

to run a named tests:
```
CoGAPS::run_catch_unit_tests("[tag1],[tag2]", [reporter])
```
where [tag1], [tag2] are he cpp tags assigned to test functions in the function
definitions as defined in (src/cpp_tests/[file].cpp)

## Rcpp::compileAttributes

Run every time if headers are changed.

## autoconf
usage
```
autoconf
./configure --enable-cpp-tests
```


