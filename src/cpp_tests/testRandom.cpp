#include "catch.h"
#include "../Random.h"

TEST_CASE("Test Random.h - Random Number Generation")
{
    Random::setSeed(0);

    SECTION("Make sure uniform01 is working")
    {
        REQUIRE(Random::uniform() != Random::uniform());
    }

    SECTION("Test uniform distribution over unit interval")
    {
        double min = 1, max = 0;
        for (unsigned i = 0; i < 10000; ++i)
        {
            min = std::min(Random::uniform(), min);
            max = std::max(Random::uniform(), max);
        }
        REQUIRE(min >= 0);
        REQUIRE(max <= 1);
    }

    SECTION("Test uniform distribution over general interval")
    {
        // invalid bounds
        REQUIRE_THROWS(Random::uniform(2.0, 1.4));

        // bounds equal
        REQUIRE(Random::uniform(4.3,4.3) == 4.3);

        // full range possible
        double min = 10., max = 0.;
        for (unsigned i = 0; i < 1000; ++i)
        {
            min = std::min(Random::uniform(0.0,10.0), min);
            max = std::max(Random::uniform(0.0,10.0), max);
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
            norm[i] = Random::normal(0, 1);        
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
            double num = Random::poisson(4);
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
            double num = Random::exponential(1);
            total += num;

            REQUIRE(num >= 0.0); // should be non-negative
        }
        double mean = total / 10000;
        REQUIRE(mean == Approx(1).epsilon(0.025));
    }
}

TEST_CASE("Test Random.h - Distribution Calculations")
{
    /*SECTION("Test p_exp")
    {
        REQUIRE(Random::p_exp(1, 3)    == Approx(0.9502).epsilon(0.001));
        REQUIRE(Random::p_exp(2, 1)    == Approx(0.8647).epsilon(0.001));
        REQUIRE(Random::p_exp(0.01, 1) == Approx(0.0100).epsilon(0.001));
        REQUIRE(Random::p_exp(0, 5)    == Approx(0.0000).epsilon(0.001));

        REQUIRE_THROWS(Random::p_exp(1, 0));
    }

    SECTION("Test q_exp")
    {
        REQUIRE(Random::q_exp(0, 2)    == Approx(0.0000).epsilon(0.001));
        REQUIRE(Random::q_exp(0.1, 1)  == Approx(0.1054).epsilon(0.001));
        REQUIRE(Random::q_exp(0.25, 1) == Approx(0.2877).epsilon(0.001));
        REQUIRE(Random::q_exp(0.5, 1)  == Approx(0.6931).epsilon(0.001));
        REQUIRE(Random::q_exp(0.75, 1) == Approx(1.3863).epsilon(0.001));
        REQUIRE(Random::q_exp(0.99, 1) == Approx(4.6052).epsilon(0.001));

        REQUIRE_THROWS(Random::q_exp(0.5, 0));
        REQUIRE_THROWS(Random::q_exp(-1, 1));
        REQUIRE_THROWS(Random::q_exp(1.1, 1));
        REQUIRE_THROWS(Random::q_exp(1, 1));
    }*/

    SECTION("Test d_gamma")
    {
        REQUIRE(Random::d_gamma(0.5, 1, 1) == Approx(0.6065).epsilon(0.001));
    }

    SECTION("Test p_gamma")
    {
        REQUIRE(Random::p_gamma(0.5, 1, 1) == Approx(0.3935).epsilon(0.001));
    }


    SECTION("Test q_gamma")
    {
        REQUIRE(Random::q_gamma(0.5, 1, 1) == Approx(0.6931).epsilon(0.001));
    }


    SECTION("Test d_norm")
    {
        REQUIRE(Random::d_norm(0.5, 0, 1) == Approx(0.3521).epsilon(0.001));
    }


    SECTION("Test q_norm")
    {
        REQUIRE(Random::q_norm(0.5, 0, 1) == Approx(0.0000).epsilon(0.001));
    }


    SECTION("Test p_norm")
    {
        REQUIRE(Random::p_norm(0.5, 0, 1) == Approx(0.6915).epsilon(0.001));
    }
}

