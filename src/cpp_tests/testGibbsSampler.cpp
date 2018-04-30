#include "catch.h"

#include "../GibbsSampler.h"
#include "../Algorithms.h"

#include <Rcpp.h>

#define TEST_APPROX(x) Approx(x).epsilon(0.001)

#if 0

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
        GibbsSampler sampler(rD, rS, 5, 0.01f, 0.01f, 1.f, 1.f, false,
            'N', Rcpp::NumericMatrix(), PUMP_UNIQUE);

        REQUIRE(sampler.totalNumAtoms('A') == 0);
        REQUIRE(sampler.totalNumAtoms('P') == 0);

        TwoWayMatrix D(rD), S(rS), AP(1363, 9);
        REQUIRE(sampler.chi2() == 2.f * gaps::algo::loglikelihood(D, S, AP));
    }

    SECTION("Create GibbsSampler with Fixed Patterns")
    {
        //TODO
    }

    SECTION("Update GibbsSampler")
    {
        GibbsSampler sampler(rD, rS, 5, 0.01f, 0.01f, 1.f, 1.f, false,
            'N', Rcpp::NumericMatrix(), PUMP_UNIQUE);

        for (unsigned i = 0; i < 1000; ++i)
        {
            REQUIRE_NOTHROW(sampler.update('A'));
            REQUIRE_NOTHROW(sampler.update('P'));
        }
        REQUIRE(sampler.totalNumAtoms('A') > 0);
        REQUIRE(sampler.totalNumAtoms('P') > 0);
    }

    SECTION("GibbsSampler Statistics")
    {
        GibbsSampler sampler(rD, rS, 5, 0.01f, 0.01f, 1.f, 1.f, false,
            'N', Rcpp::NumericMatrix(), PUMP_UNIQUE);

        for (unsigned i = 0; i < 1000; ++i)
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
        REQUIRE(AMean.ncol() == 5);

        REQUIRE(AStd.nrow() == rD.nrow());
        REQUIRE(AStd.ncol() == 5);

        REQUIRE(PMean.nrow() == 5);
        REQUIRE(PMean.ncol() == rD.ncol());

        REQUIRE(PStd.nrow() == 5);
        REQUIRE(PStd.ncol() == rD.ncol());

        for (signed r = 0; r < AMean.nrow(); ++r)
        {
            for (signed c = 0; c < AMean.ncol(); ++c)
            {
                REQUIRE(AMean(r,c) >= 0.0);
                REQUIRE(AStd(r,c) >= 0.0);
            }
        }

        for (signed r = 0; r < PMean.nrow(); ++r)
        {
            for (signed c = 0; c < PMean.ncol(); ++c)
            {
                REQUIRE(PMean(r,c) >= 0.0);
                REQUIRE(PStd(r,c) >= 0.0);
            }
        }
        REQUIRE(sampler.meanChiSq() == TEST_APPROX(8371.568));
    }
}

#ifdef GAPS_INTERNAL_TESTS

TEST_CASE("Internal GibbsSampler Tests")
{
    gaps::random::setSeed(0);

    Rcpp::Function asMatrix("as.matrix");
    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CoGAPS");
    Rcpp::NumericMatrix rD = asMatrix(pkgEnv.find("GIST.D"));
    Rcpp::NumericMatrix rS = asMatrix(pkgEnv.find("GIST.D"));

    SECTION("Test deltaLL and chi2 consistency")
    {
        GibbsSampler sampler(rD, rS, 5, 0.01f, 0.01f, 1.f, 1.f, false,
            'N', Rcpp::NumericMatrix());

        // need to populate the matrices a little bit
        for (unsigned i = 0; i < 5000; ++i)
        {
            REQUIRE_NOTHROW(sampler.update('A'));
            REQUIRE_NOTHROW(sampler.update('P'));
        }

        for (unsigned i = 0; i < 5000; ++i)
        {
            float preLL = sampler.chi2();
            AtomicProposal prop = sampler.mADomain.makeProposal();
            MatrixChange ch = sampler.mADomain.acceptProposal(prop);
            float delLL = sampler.computeDeltaLL(ch);
            if (std::abs(delLL) < 1E6) // large values == large truncation error
            {
                sampler.mAMatrix.update(ch);
                sampler.updateAPMatrix(ch);
                REQUIRE(preLL - 2.f * delLL == TEST_APPROX(sampler.chi2()));
            }
        }

        for (unsigned i = 0; i < 5000; ++i)
        {
            float preLL = sampler.chi2();
            AtomicProposal prop = sampler.mPDomain.makeProposal();
            MatrixChange ch = sampler.mPDomain.acceptProposal(prop);
            float delLL = sampler.computeDeltaLL(ch);
            if (std::abs(delLL) < 1E6) // large values == large truncation error
            {
                sampler.mPMatrix.update(ch);
                sampler.updateAPMatrix(ch);
                REQUIRE(preLL - 2.f * delLL == TEST_APPROX(sampler.chi2()));
            }
        }
    }
}

#endif

#endif