#include "catch.h"
#include "../Random.h"

TEST_CASE("Test Random.h")
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

