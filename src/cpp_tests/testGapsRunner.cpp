#include "catch.h"
#include "../GapsRunner.h"

TEST_CASE("Test Top Level CoGAPS Call")
{
    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CoGAPS");
    Rcpp::NumericMatrix D = pkgEnv.find("SimpSim.D");
    Rcpp::NumericMatrix S = pkgEnv.find("SimpSim.S");

    GapsRunner runner(D, S, 3, 500, 100, 500, 250, 0, 0.01f, 0.01f, 100.f,
        100.f, 123, false, false, 0, "gaps_checkpoint.out", 'N',
        Rcpp::NumericMatrix(1,1));
    REQUIRE_NOTHROW(runner.run());
}

#ifdef GAPS_INTERNAL_TESTS

TEST_CASE("Test Individual Steps")
{

}

#endif