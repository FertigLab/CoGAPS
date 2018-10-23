#include "catch.h"
#include "../data_structures/Vector.h"
#include "../math/Random.h"
#include "../math/VectorMath.h"

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