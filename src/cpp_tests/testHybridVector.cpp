#include "catch.h"
#include "../data_structures/HybridVector.h"
#include "../math/Random.h"
#include "../math/Math.h"
#include "../math/VectorMath.h"

TEST_CASE("Test HybridVector.h")
{
    GapsRandomState randState(123);

    SECTION("Test size constructor")
    {
        HybridVector v(100);
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
        HybridVector v(in_v);

        REQUIRE(v.size() == 1000);
        REQUIRE(!gaps::isVectorZero(v));
        REQUIRE(gaps::max(v) <= 1.f);
        REQUIRE(gaps::min(v) >= 0.f);
    }

    SECTION("TEST add function")
    {
        // create vector
        GapsRng rng(&randState);
        HybridVector v(1000);

        // add tons of values to it, occasionally erasing
        for (unsigned n = 0; n < 100000; ++n)
        {
            unsigned ndx = rng.uniform32(0,999);
            float diff = rng.uniform();
            float old = v[ndx];
            if (rng.uniform() < 0.2f)
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

    #ifdef GAPS_INTERNAL_TESTS
        for (unsigned i = 0; i < v.size(); ++i)
        {
            if (v[i] == 0.f)
            {
                uint64_t mask = v.mIndexBitFlags[i / 64] & (1ull << (i % 64));
                REQUIRE(!mask);
            }
            else
            {
                uint64_t mask = v.mIndexBitFlags[i / 64] & (1ull << (i % 64));
                REQUIRE(mask);
            }
        }
    #endif
    }
}