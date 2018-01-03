#include "catch.h"

#include "../GibbsSampler.h"

#include <Rcpp.h>

TEST_CASE("Test GibbsSampler.h")
{
/*    gaps::random::setSeed(0);

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
        GibbsSampler sampler(rD, rS, 10, 0.01, 0.01, 1.0, 1.0, false);

        REQUIRE(sampler.chi2() == 24534.0);
        REQUIRE(sampler.totalNumAtoms('A') == 0);
        REQUIRE(sampler.totalNumAtoms('P') == 0);
    }

    SECTION("Update GibbsSampler")
    {
        GibbsSampler sampler(rD, rS, 10, 0.01, 0.01, 1.0, 1.0, false);

        for (unsigned i = 0; i < 100; ++i)
        {
            REQUIRE_NOTHROW(sampler.update('A'));
            REQUIRE_NOTHROW(sampler.update('P'));
        }
        REQUIRE(sampler.totalNumAtoms('A') > 10);
        REQUIRE(sampler.totalNumAtoms('P') > 10);
    }

    SECTION("GibbsSampler Statistics")
    {
        GibbsSampler sampler(rD, rS, 10, 0.01, 0.01, 1.0, 1.0, false);

        for (unsigned i = 0; i < 100; ++i)
        {
            for (unsigned j = 0; j < 10; ++j)
            {
                sampler.update('A');
                sampler.update('P');
            }
            sampler.updateStatistics();
        }
        Rcpp::NumericMatrix AMean = sampler.AMeanRMatrix();
        Rcpp::NumericMatrix AStd = sampler.AStdRMatrix();
        Rcpp::NumericMatrix PMean = sampler.PMeanRMatrix();
        Rcpp::NumericMatrix PStd = sampler.PStdRMatrix();

        REQUIRE(AMean.nrow() == rD.nrow());
        REQUIRE(AMean.ncol() == 10);

        REQUIRE(AStd.nrow() == rD.nrow());
        REQUIRE(AStd.ncol() == 10);

        REQUIRE(PMean.nrow() == 10);
        REQUIRE(PMean.ncol() == rD.ncol());

        REQUIRE(PStd.nrow() == 10);
        REQUIRE(PStd.ncol() == rD.ncol());

        for (unsigned r = 0; r < AMean.nrow(); ++r)
        {
            for (unsigned c = 0; c < AMean.ncol(); ++c)
            {
                REQUIRE(AMean(r,c) >= 0.0);
                REQUIRE(AStd(r,c) >= 0.0);
            }
        }

        for (unsigned r = 0; r < PMean.nrow(); ++r)
        {
            for (unsigned c = 0; c < PMean.ncol(); ++c)
            {
                REQUIRE(PMean(r,c) >= 0.0);
                REQUIRE(PStd(r,c) >= 0.0);
            }
        }
    }
*/
}

#ifdef GAPS_INTERNAL_TESTS

TEST_CASE("Internal GibbsSampler Tests")
{

}

#endif