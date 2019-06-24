#include "catch.h"
#include "../data_structures/SparseVector.h"
#include "../math/Random.h"
#include "../math/Math.h"
#include "../math/VectorMath.h"

TEST_CASE("Test SparseVector.h")
{
    GapsRandomState randState(123);

    SECTION("Test size constructor")
    {
        SparseVector sv(100);
        REQUIRE(sv.size() == 100);

        Vector v(sv.getDense());
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
        SparseVector sv(in_v);
        REQUIRE(sv.size() == 1000);

        Vector v(sv.getDense());
        REQUIRE(v.size() == 1000);
        REQUIRE(!gaps::isVectorZero(v));
        REQUIRE(gaps::max(v) <= 1.f);
        REQUIRE(gaps::min(v) >= 0.f);
    }

    SECTION("Test Vector constructor")
    {
        GapsRng rng(&randState);
        std::vector<float> in_v;
        for (unsigned n = 0; n < 1000; ++n)
        {
            in_v.push_back(rng.uniform());
        }

        Vector v1(in_v);
        SparseVector sv(v1);
        Vector v2(sv.getDense());

        REQUIRE(v1.size() == sv.size());
        REQUIRE(sv.size() == v2.size());

        for (unsigned i = 0; i < v1.size(); ++i)
        {
            REQUIRE(v1[i] == v2[i]);
        }
    }

#if 0
    SECTION("bit flags set correctly")
    {
        SparseVector v(10);
        v.insert(0, 1.f);
        v.insert(4, 5.f);
        v.insert(7, 8.f);
        v.insert(9, 10.f);
        REQUIRE(v.mIndexBitFlags[0] == 0b1010010001);
    }

    SECTION("values placed correctly")
    {
        SparseVector v(10);
        v.insert(0, 1.f);
        v.insert(4, 5.f);
        v.insert(7, 8.f);
        v.insert(9, 10.f);

        REQUIRE(v.mData.size() == 4);
        REQUIRE(v.mData[0] == 1.f);
        REQUIRE(v.mData[1] == 5.f);
        REQUIRE(v.mData[2] == 8.f);
        REQUIRE(v.mData[3] == 10.f);
    }
#endif
}