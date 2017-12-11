#include "catch.h"

#include "../GibbsSampler.h"

#include <Rcpp.h>

TEST_CASE("Test GibbsSampler.h")
{
    gaps::random::setSeed(0);

    Rcpp::Function asMatrix("as.matrix");
    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CoGAPS");
    Rcpp::NumericMatrix rD = asMatrix(pkgEnv.find("GIST.D"));
    Rcpp::NumericMatrix rS = asMatrix(pkgEnv.find("GIST.D"));

    REQUIRE(rD.nrow() == 1363);
    REQUIRE(rD.ncol() == 9);

    REQUIRE(rS.nrow() == 1363);
    REQUIRE(rS.ncol() == 9);

    SECTION("Create GibbsSampler")
    {
        GibbsSampler sampler(rD, rS, 10, 1, 1, 5.0, 5.0, false);

        REQUIRE(sampler.chi2() == 0.0);
        REQUIRE(sampler.totalNumAtoms('A') == 0);
        REQUIRE(sampler.totalNumAtoms('P') == 0);
    }

    SECTION("Update GibbsSampler")
    {
        GibbsSampler sampler(rD, rS, 10, 1, 1, 5.0, 5.0, false);

        for (unsigned i = 0; i < 1000; ++i)
        {
            REQUIRE_NOTHROW(sampler.update('A'));
            REQUIRE_NOTHROW(sampler.update('P'));
        }
    }
}