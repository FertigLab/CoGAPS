#include <testthat.h>
#include "../testthat-tweak.h"
#include "../math/Random.h"
#include "../math/Math.h"
#include "../utils/GapsPrint.h"

#include <algorithm>
#include <cmath>

#define TEST_APPROX(x) Approx(x).epsilon(0.001)

// this is intended to replicated the random stream that happens when
// each proposal is creating a new rng and using it a few times

GapsRandomState testRandState(42);
GapsRng testRng(&testRandState);

TEST_CASE("Random Number Generation -- basic","[randomrng][basic]")
{
    SECTION("Make sure uniform is working")
    {
        REQUIRE(testRng.uniform64() != testRng.uniform64());
        REQUIRE(testRng.uniform64(0,1000) != testRng.uniform64(0,1000));
        REQUIRE(testRng.uniform64(1000,1000) == 1000);
        REQUIRE(testRng.uniform64(0,2)<3);
        REQUIRE(testRng.uniform32() != testRng.uniform32());
        REQUIRE(testRng.uniform32(0,100) != testRng.uniform32(0,100));
        REQUIRE(testRng.uniform32(1000,1000) == 1000);
    }
    SECTION("Make sure uniform is working on size_t")
    {
        size_t zero=0, thou=1000;
        REQUIRE(testRng.uniform64(zero,thou) != testRng.uniform64(zero,thou));
        REQUIRE(testRng.uniform64(thou,thou) == thou);
    }
}

static void requireSmallError(float in, float out, float est, float tol)
{
    float denom = std::max(std::abs(out), 1.f);
    if (std::abs(est - out) / denom >= tol)
    {
        gaps_printf("input: %f, output: %f, error: %f\n", in, out,
            std::abs(est - out));
    }
    REQUIRE(std::abs(est - out) / denom < tol);
}

TEST_CASE("Test error of q_norm lookup table","[randomrng][q_norm]")
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

TEST_CASE("Test error of p_norm lookup table","[randomrng][p_norm]")
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

// ported from the removed gaps::random:: global API to the current GapsRng object.
// GapsRng is deterministic for a fixed seed, so these are stable, not flaky.
TEST_CASE("Random.h - RNG distributions","[randomrng][distributions]")
{
    GapsRandomState randState(0);
    GapsRng rng(&randState);

    SECTION("uniform01 produces varying values")
    {
        REQUIRE(rng.uniform() != rng.uniform());
    }

    SECTION("uniform over the unit interval")
    {
        float mn = 1.f, mx = 0.f, sum = 0.f;
        const unsigned N = 10000;
        for (unsigned i = 0; i < N; ++i)
        {
            float u = rng.uniform();
            mn = u < mn ? u : mn;
            mx = u > mx ? u : mx;
            sum += u;
        }
        REQUIRE(sum / N == Approx(0.5f).epsilon(0.02f));
        REQUIRE(mn >= 0.f);
        REQUIRE(mn < 0.02f);
        REQUIRE(mx <= 1.f);
        REQUIRE(mx > 0.98f);
    }

    SECTION("uniform over a general interval")
    {
        REQUIRE(rng.uniform(4.3f, 4.3f) == 4.3f);
        float mn = 10.f, mx = 0.f;
        for (unsigned i = 0; i < 1000; ++i)
        {
            float u = rng.uniform(0.f, 10.f);
            mn = u < mn ? u : mn;
            mx = u > mx ? u : mx;
        }
        REQUIRE(mn < 0.2f);
        REQUIRE(mx > 9.8f);
    }

    SECTION("poisson mean")
    {
        double total = 0.;
        const unsigned N = 10000;
        for (unsigned i = 0; i < N; ++i)
        {
            int num = rng.poisson(4.0);
            REQUIRE(num >= 0);
            total += num;
        }
        REQUIRE(total / N == Approx(4.0).epsilon(0.03));
    }

    SECTION("exponential mean")
    {
        double total = 0.;
        const unsigned N = 10000;
        for (unsigned i = 0; i < N; ++i)
        {
            float num = rng.exponential(1.f);
            REQUIRE(num >= 0.f);
            total += num;
        }
        REQUIRE(total / N == Approx(1.0).epsilon(0.03));
    }
}

TEST_CASE("Random.h - distribution calculations","[randomrng][distcalc]")
{
    REQUIRE(gaps::d_gamma(0.5f, 1.f, 1.f) == TEST_APPROX(0.607f));
    REQUIRE(gaps::p_gamma(0.5f, 1.f, 1.f) == TEST_APPROX(0.394f));
    REQUIRE(gaps::q_gamma(0.5f, 1.f, 1.f) == TEST_APPROX(0.693f));
    REQUIRE(gaps::d_norm(0.5f, 0.f, 1.f) == TEST_APPROX(0.352f));
    REQUIRE(gaps::q_norm(0.5f, 0.f, 1.f) == TEST_APPROX(0.000f));
    REQUIRE(gaps::p_norm(0.5f, 0.f, 1.f) == TEST_APPROX(0.692f));
}

// Tooling harness (not a unit test): writes a long random stream to a file for
// external diehard RNG tests. Left disabled; enable manually when needed.
#if 0
#include "../utils/Archive.h"

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
