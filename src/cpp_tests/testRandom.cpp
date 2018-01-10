#include "catch.h"
#include "../Random.h"

TEST_CASE("Test Random.h - Random Number Generation")
{
    gaps::random::setSeed(0);

    SECTION("Make sure uniform01 is working")
    {
        REQUIRE(gaps::random::uniform() != gaps::random::uniform());
    }

    SECTION("Test uniform distribution over unit interval")
    {
        double min = 1, max = 0;
        double sum = 0.0;
        unsigned N = 10000;
        for (unsigned i = 0; i < N; ++i)
        {
            min = std::min(gaps::random::uniform(), min);
            max = std::max(gaps::random::uniform(), max);
            sum += gaps::random::uniform();
        }
        REQUIRE(sum / N == Approx(0.5).epsilon(0.01));
        REQUIRE(min >= 0);
        REQUIRE(min < 0.01);
        REQUIRE(max <= 1);
        REQUIRE(max > 0.99);
    }

    SECTION("Test uniform distribution over general interval")
    {
        // bounds equal
        REQUIRE(gaps::random::uniform(4.3,4.3) == 4.3);

        // full range possible
        double min = 10., max = 0.;
        for (unsigned i = 0; i < 1000; ++i)
        {
            min = std::min(gaps::random::uniform(0.0,10.0), min);
            max = std::max(gaps::random::uniform(0.0,10.0), max);
        }
        REQUIRE(min < 0.1);
        REQUIRE(max > 9.9);
    }

    SECTION("Test uniform distribution over integer range")
    {
        // TODO
    }

    SECTION("Test uniform distribution over 64 bit integers")
    {
        // TODO
    }

    SECTION("Test normal distribution")
    {
        // sample distribution
        double mean = 0.0, var = 0.0;
        double norm[1024];
        for (unsigned i = 0; i < 1024; ++i)
        {
            norm[i] = gaps::random::normal(0, 1);        
            mean += norm[i];
        }

        // check parameters
        mean /= 1024;
        for (unsigned i = 0; i < 1000; ++i)
        {
            var += pow(norm[i] - mean, 2);
        }
        var /= 1024;
        REQUIRE(mean == Approx(0).epsilon(0.025));
        REQUIRE(var == Approx(1).epsilon(0.025));
    }

    SECTION("Test poisson distribution")
    {
        double total = 0;
        for (unsigned i = 0; i < 10000; ++i)
        {
            double num = gaps::random::poisson(4);
            total += num;

            REQUIRE((int)num == num); // should be integer
            REQUIRE(num >= 0.0); // should be non-negative
        }
        double mean = total / 10000;
        REQUIRE(mean == Approx(4).epsilon(0.025));
    }

    SECTION("Test exponential distribution")
    {
        double total = 0;
        for (unsigned i = 0; i < 10000; ++i)
        {
            double num = gaps::random::exponential(1);
            total += num;

            REQUIRE(num >= 0.0); // should be non-negative
        }
        double mean = total / 10000;
        REQUIRE(mean == Approx(1).epsilon(0.025));
    }
}

TEST_CASE("Test Random.h - Distribution Calculations")
{
    SECTION("Test d_gamma")
    {
        REQUIRE(gaps::random::d_gamma(0.5, 1, 1) == Approx(0.6065).epsilon(0.001));
    }

    SECTION("Test p_gamma")
    {
        REQUIRE(gaps::random::p_gamma(0.5, 1, 1) == Approx(0.3935).epsilon(0.001));
    }


    SECTION("Test q_gamma")
    {
        REQUIRE(gaps::random::q_gamma(0.5, 1, 1) == Approx(0.6931).epsilon(0.001));
    }


    SECTION("Test d_norm")
    {
        REQUIRE(gaps::random::d_norm(0.5, 0, 1) == Approx(0.3521).epsilon(0.001));
    }


    SECTION("Test q_norm")
    {
        REQUIRE(gaps::random::q_norm(0.5, 0, 1) == Approx(0.0000).epsilon(0.001));
    }


    SECTION("Test p_norm")
    {
        REQUIRE(gaps::random::p_norm(0.5, 0, 1) == Approx(0.6915).epsilon(0.001));
    }
}

