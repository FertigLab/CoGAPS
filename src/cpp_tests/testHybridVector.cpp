#include "catch.h"
#include "../data_structures/HybridVector.h"
#include "../math/Random.h"
#include "../math/Math.h"
#include "../math/VectorMath.h"

TEST_CASE("Test HybridVector.h")
{
    GapsRng::setSeed(123);

    SECTION("Test size constructor")
    {
        HybridVector v(100);
        REQUIRE(v.size() == 100);
        REQUIRE(gaps::isVectorZero(v));
        REQUIRE(gaps::sum(v) == 0.f);
    }

    SECTION("Test std::vector constructor")
    {   
        GapsRng rng;
        std::vector<float> in_v;
        for (unsigned n = 0; n < 1000; ++n)
        {
            in_v.push_back(rng.uniform());
        }
        HybridVector v(in_v);

        REQUIRE(v.size() == 1000);
        REQUIRE(!gaps::isVectorZero(v));
        REQUIRE(gaps::max(v) <= 1.f);
        REQUIRE(gaps::min(v) >= 0.f);
    }

    SECTION("TEST add function")
    {
        // create vector
        GapsRng rng;
        HybridVector v(1000);

        // add tons of values to it, occasionally erasing
        for (unsigned n = 0; n < 100000; ++n)
        {
            unsigned ndx = rng.uniform32(0,999);
            float diff = rng.uniform();
            float old = v[ndx];
            if (rng.uniform() < 0.05f)
            {
                diff = -1.f * old;
                REQUIRE(v.add(ndx, diff));
            }
            else
            {
                REQUIRE(!v.add(ndx, diff));
                REQUIRE(v[ndx] == old + diff);
            }
        }
    }
}