#include "catch.h"
#include "../data_structures/HashSets.h"
#include "../math/Random.h"

TEST_CASE("Test HashSets.h - FixedHashSetU32")
{
    GapsRandomState randState(123);

    FixedHashSetU32 hSet(1000);
    REQUIRE(hSet.isEmpty());

    // burn in
    GapsRng rng(&randState);
    for (unsigned n = 0; n < 1000; ++n)
    {
        uint64_t u = 0;
        for (unsigned i = 0; i < 100; ++i)
        {
            u = rng.uniform32(0, 1000);
            hSet.insert(u);
            REQUIRE(hSet.contains(u));
        }
        hSet.clear();
        REQUIRE(!hSet.contains(u));
        REQUIRE(hSet.isEmpty());
    }
}

TEST_CASE("Test HashSets.h - SmallHashSetU64")
{
    GapsRandomState randState(123);

    SmallHashSetU64 hSet;
    REQUIRE(hSet.isEmpty());

    // burn in
    GapsRng rng(&randState);
    for (unsigned n = 0; n < 1000; ++n)
    {
        uint64_t u = 0;
        for (unsigned i = 0; i < 100; ++i)
        {
            u = rng.uniform64();
            hSet.insert(u);
            REQUIRE(hSet.contains(u));
        }
        hSet.clear();
        REQUIRE(!hSet.contains(u));
        REQUIRE(hSet.isEmpty());
    }
}

TEST_CASE("Test HashSets.h - SmallPairedHashSetU64")
{
    SmallPairedHashSetU64 hSet;
    REQUIRE(hSet.isEmpty());

    GapsRandomState randState(123);
    GapsRng rng(&randState);
    for (unsigned n = 0; n < 1000; ++n)
    {
        for (unsigned i = 0; i < 100; ++i)
        {
            uint64_t u1 = rng.uniform64();
            uint64_t u2 = rng.uniform64();
            hSet.insert(u1, u2);
            REQUIRE(hSet.contains(u1));
            REQUIRE(hSet.contains(u2));
            uint64_t d = u1 > u2 ? u1 - u2 : u2 - u1;
            REQUIRE(hSet.overlap(u1 > u2 ? u2 + d/2 : u1 + d/2));
        }
        hSet.clear();
        REQUIRE(hSet.isEmpty());
    }
}


