#include "catch.h"
#include "../randgen.h"
#include "../sub_func.h"

TEST_CASE("Test randgen.h")
{
    rng.seed(0);

    SECTION("Test uniform distribution over unit interval")
    {
        double min = 1, max = 0;
        for (unsigned i = 0; i < 10000; ++i)
        {
            min = std::min(randgen('U'), min);
            max = std::max(randgen('U'), max);
        }
        REQUIRE(min >= 0);
        REQUIRE(max <= 1);
    }

    SECTION("Test uniform distribution over general interval")
    {
        // invalid bounds
        REQUIRE_THROWS(gaps::sub_func::runif(2.0, 1.4));

        // bounds equal
        REQUIRE(gaps::sub_func::runif(4.3,4.3) == 4.3);

        // full range possible
        double min = 10., max = 0.;
        for (unsigned i = 0; i < 1000; ++i)
        {
            min = std::min(gaps::sub_func::runif(0.0,10.0), min);
            max = std::max(gaps::sub_func::runif(0.0,10.0), max);
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
            norm[i] = randgen('N', 0, 1);        
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
            double num = randgen('P', 4);
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
            double num = randgen('E', 1);
            total += num;

            REQUIRE(num >= 0.0); // should be non-negative
        }
        double mean = total / 10000;
        REQUIRE(mean == Approx(1).epsilon(0.025));
    }
}

