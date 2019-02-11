#include "catch.h"
#include "../data_structures/Vector.h"
#include "../math/Random.h"
#include "../math/VectorMath.h"

// optional test used for benchmarking, set to 0 to disable, 1 to enable
#if 0

// boost time helpers
#include <boost/date_time/posix_time/posix_time.hpp>
namespace bpt = boost::posix_time;
#define bpt_now() bpt::microsec_clock::local_time()

TEST_CASE("Benchmark Dot Product")
{
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    std::vector<Vector> mVecs;
    for (unsigned i = 0; i < 1300; ++i)
    {
        mVecs.push_back(Vector(50000));
        for (unsigned j = 0; j < mVecs[i].size(); ++j)
        {
            mVecs[i][j] = rng.uniform(0.f, 100.f);
        }
    }

    float sum = 0.f;
    bpt::ptime start = bpt_now();
    for (unsigned i = 0; i < mVecs.size(); ++i)
    {
        for (unsigned j = i; j < mVecs.size(); ++j)
        {
            sum += gaps::dot(mVecs[i], mVecs[j]);
        }
    }
    bpt::time_duration diff = bpt_now() - start;
    gaps_printf("-------\n-------\n-------\n-------\n", sum);
    gaps_printf("sum: %f\n", sum);
    gaps_printf("dot product milliseconds: %lu\n", diff.total_milliseconds());
    gaps_printf("-------\n-------\n-------\n-------\n", sum);
}
#endif

TEST_CASE("Test Vector.h")
{
    GapsRandomState randState(123);

    SECTION("Test size constructor")
    {
        Vector v(100);
        REQUIRE(v.size() == 100);
        REQUIRE(gaps::isVectorZero(v));
        REQUIRE(gaps::sum(v) == 0.f);
    }

    SECTION("Test std::vector constructor")
    {   
        GapsRng rng(&randState);
        std::vector<float> in_v;
        for (unsigned n = 0; n < 1000; ++n)
        {
            in_v.push_back(rng.uniform());
        }
        Vector v(in_v);

        REQUIRE(v.size() == 1000);
        REQUIRE(!gaps::isVectorZero(v));
        REQUIRE(gaps::max(v) <= 1.f);
        REQUIRE(gaps::min(v) >= 0.f);
    }

    SECTION("TEST += operator")
    {
        GapsRng rng(&randState);
        std::vector<float> in_v;
        for (unsigned n = 0; n < 1000; ++n)
        {
            in_v.push_back(rng.uniform());
        }
        Vector v(in_v);

        float s = gaps::sum(v);
        v += v;
        REQUIRE(gaps::sum(v) == 2.f * s);
    }
}