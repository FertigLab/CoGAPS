#include "catch.h"
#include "../math/Random.h"

#define TEST_APPROX(x) Approx(x).epsilon(0.001)

TEST_CASE("Test Random.h - Random Number Generation")
{
    gaps::random::Generator::setSeed(0);

    SECTION("Make sure uniform01 is working")
    {
        REQUIRE(gaps::random::Generator::uniform() != gaps::random::Generator::uniform());
    }

    SECTION("Test uniform distribution over unit interval")
    {
        float min = 1.f, max = 0.f;
        float sum = 0.f;
        unsigned N = 10000;
        for (unsigned i = 0; i < N; ++i)
        {
            min = std::min(gaps::random::Generator::uniform(), min);
            max = std::max(gaps::random::Generator::uniform(), max);
            sum += gaps::random::Generator::uniform();
        }
        REQUIRE(sum / N == Approx(0.5f).epsilon(0.01f));
        REQUIRE(min >= 0.f);
        REQUIRE(min < 0.01f);
        REQUIRE(max <= 1.f);
        REQUIRE(max > 0.99f);
    }

    SECTION("Test uniform distribution over general interval")
    {
        // bounds equal
        REQUIRE(gaps::random::Generator::uniform(4.3f, 4.3f) == 4.3f);

        // full range possible
        float min = 10., max = 0.;
        for (unsigned i = 0; i < 1000; ++i)
        {
            min = std::min(gaps::random::Generator::uniform(0.f,10.f), min);
            max = std::max(gaps::random::Generator::uniform(0.f,10.f), max);
        }
        REQUIRE(min < 0.1f);
        REQUIRE(max > 9.9f);
    }

    SECTION("Test uniform distribution over 64 bit integers")
    {
        // TODO
    }

    SECTION("Test normal distribution")
    {
        // sample distribution
        float mean = 0.f, var = 0.f;
        float norm[1024];
        for (unsigned i = 0; i < 1024; ++i)
        {
            norm[i] = gaps::random::Generator::normal(0.f, 1.f);
            mean += norm[i];
        }

        // check parameters
        mean /= 1024.f;
        for (unsigned i = 0; i < 1000; ++i)
        {
            var += pow(norm[i] - mean, 2);
        }
        var /= 1024.f;
        REQUIRE(mean == Approx(0.f).epsilon(0.025f));
        REQUIRE(var == Approx(1.f).epsilon(0.025f));
    }

    SECTION("Test poisson distribution")
    {
        float total = 0.f;
        for (unsigned i = 0; i < 10000; ++i)
        {
            float num = gaps::random::Generator::poisson(4.f);
            total += num;

            REQUIRE((int)num == num); // should be integer
            REQUIRE(num >= 0.f); // should be non-negative
        }
        float mean = total / 10000.f;
        REQUIRE(mean == Approx(4.f).epsilon(0.025f));
    }

    SECTION("Test exponential distribution")
    {
        float total = 0.f;
        for (unsigned i = 0; i < 10000; ++i)
        {
            float num = gaps::random::Generator::exponential(1.f);
            total += num;

            REQUIRE(num >= 0.f); // should be non-negative
        }
        float mean = total / 10000.f;
        REQUIRE(mean == Approx(1.f).epsilon(0.025f));
    }
}

TEST_CASE("Test Random.h - Distribution Calculations")
{
    REQUIRE(gaps::random::d_gamma(0.5f, 1.f, 1.f) == TEST_APPROX(0.607f));
    REQUIRE(gaps::random::p_gamma(0.5f, 1.f, 1.f) == TEST_APPROX(0.394f));
    REQUIRE(gaps::random::q_gamma(0.5f, 1.f, 1.f) == TEST_APPROX(0.693f));
    REQUIRE(gaps::random::d_norm(0.5f, 0.f, 1.f) == TEST_APPROX(0.352f));
    REQUIRE(gaps::random::q_norm(0.5f, 0.f, 1.f) == TEST_APPROX(0.000f));
    REQUIRE(gaps::random::p_norm(0.5f, 0.f, 1.f) == TEST_APPROX(0.692f));
}
