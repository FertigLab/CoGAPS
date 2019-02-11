#include "catch.h"
#include "../math/Random.h"
#include "../math/Math.h"

#define TEST_APPROX(x) Approx(x).epsilon(0.001)

// this is intended to replicated the random stream that happens when
// each proposal is creating a new rng and using it a few times
class EmulatedRng
{
public:

    EmulatedRng(unsigned seed)
        : randState(seed), rng(&randState), tickRng(&randState),
        remaining(tickRng.uniform32(1, 5))
    {}

    uint64_t uniform64()
    {
        advance();
        return rng.uniform64();
    }

private:

    void advance()
    {
        --remaining;
        if (remaining == 0)
        {
            rng = GapsRng(&randState);
            remaining = tickRng.uniform32(1, 5);
        }
    }

    GapsRandomState randState;
    GapsRng rng;
    GapsRng tickRng;
    unsigned remaining;
};

static void requireSmallError(float in, float out, float est, float tol)
{
    float denom = gaps::max(std::abs(out), 1.f);
    if (std::abs(est - out) / denom >= tol)
    {
        gaps_printf("input: %f, output: %f, error: %f\n", in, out,
            std::abs(est - out));
    }
    REQUIRE(std::abs(est - out) / denom < tol);
}

TEST_CASE("Test error of q_norm lookup table")
{
    GapsRandomState randState(123);

    const unsigned nIterations = 10000;
    const float mean = 0.f;
    const float sd = 1.f;
    const float tolerance = 0.03f;
    for (unsigned i = 1; i < nIterations; ++i)
    {
        float q = static_cast<float>(i) / static_cast<float>(nIterations);
        float lookup_val = randState.q_norm_fast(q, mean, sd);
        float actual_val = gaps::q_norm(q, mean, sd);
        requireSmallError(q, actual_val, lookup_val, tolerance);
    }
}

TEST_CASE("Test error of p_norm lookup table")
{
    GapsRandomState randState(123);

    const unsigned nIterations = 100 * 100; // needs to be multiple of 100
    const float mean = 0.f;
    const float sd = 1.f;
    const float tolerance = 0.03f;
    for (unsigned i = 1; i < nIterations; ++i)
    {
        float p = static_cast<float>(i) / static_cast<float>(nIterations / 100);
        float lookup_val = randState.p_norm_fast(p, mean, sd);
        float actual_val = gaps::p_norm(p, mean, sd);
        requireSmallError(p, actual_val, lookup_val, tolerance);
    }
}

#if 0
#ifdef GAPS_INTERNAL_TESTS
TEST_CASE("write random file to use in diehard tests")
{
    Archive ar("random_stream.out", ARCHIVE_WRITE);
    
    EmulatedRng rng(123);
    for (unsigned i = 0; i < 1500000; ++i)
    {
        ar << rng.uniform64();
    }
}
#endif
#endif

#if 0

TEST_CASE("Test Random.h - Random Number Generation")
{
    gaps::random::setSeed(0);

    SECTION("Make sure uniform01 is working")
    {
        REQUIRE(gaps::random::uniform() != gaps::random::uniform());
    }

    SECTION("Test uniform distribution over unit interval")
    {
        float min = 1.f, max = 0.f;
        float sum = 0.f;
        unsigned N = 10000;
        for (unsigned i = 0; i < N; ++i)
        {
            min = gaps::min(gaps::random::uniform(), min);
            max = gaps::max(gaps::random::uniform(), max);
            sum += gaps::random::uniform();
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
        REQUIRE(gaps::random::uniform(4.3f, 4.3f) == 4.3f);

        // full range possible
        float min = 10., max = 0.;
        for (unsigned i = 0; i < 1000; ++i)
        {
            min = gaps::min(gaps::random::uniform(0.f,10.f), min);
            max = gaps::max(gaps::random::uniform(0.f,10.f), max);
        }
        REQUIRE(min < 0.1f);
        REQUIRE(max > 9.9f);
    }

    SECTION("Test uniform distribution over 64 bit integers")
    {
        // TODO
    }

    SECTION("Test poisson distribution")
    {
        float total = 0.f;
        for (unsigned i = 0; i < 10000; ++i)
        {
            float num = gaps::random::poisson(4.f);
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
            float num = gaps::random::exponential(1.f);
            total += num;

            REQUIRE(num >= 0.f); // should be non-negative
        }
        float mean = total / 10000.f;
        REQUIRE(mean == Approx(1.f).epsilon(0.025f));
    }
}

TEST_CASE("Test Random.h - Distribution Calculations")
{
    REQUIRE(gaps::d_gamma(0.5f, 1.f, 1.f) == TEST_APPROX(0.607f));
    REQUIRE(gaps::p_gamma(0.5f, 1.f, 1.f) == TEST_APPROX(0.394f));
    REQUIRE(gaps::q_gamma(0.5f, 1.f, 1.f) == TEST_APPROX(0.693f));
    REQUIRE(gaps::d_norm(0.5f, 0.f, 1.f) == TEST_APPROX(0.352f));
    REQUIRE(gaps::q_norm(0.5f, 0.f, 1.f) == TEST_APPROX(0.000f));
    REQUIRE(gaps::p_norm(0.5f, 0.f, 1.f) == TEST_APPROX(0.692f));
}

#endif