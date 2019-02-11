#ifdef GAPS_CPP_UNIT_TESTS
    #define CATCH_CONFIG_RUNNER
    #include "cpp_tests/catch.h"
#endif

// [[Rcpp::export]]
int run_catch_unit_tests()
{
#ifdef GAPS_CPP_UNIT_TESTS
    Catch::Session session;
    int numFailed = session.run();
    return (numFailed < 0xFF ? numFailed : 0xFF);
#endif
    return 1;
}