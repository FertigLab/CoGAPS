#include "catch.h"
#include "../data_structures/EfficientSets.h"

TEST_CASE("Test IntFixedHashSet")
{
    IntFixedHashSet set;

    set.setDimensionSize(1000);
    REQUIRE(!set.count(123));
    
    set.insert(123);
    REQUIRE(set.count(123));

    set.clear();
    REQUIRE(!set.count(123));

    set.insert(123);
    REQUIRE(set.count(123));
}

TEST_CASE("Test IntDenseOrderedSet")
{
    IntDenseOrderedSet set;

    REQUIRE(set.isEmptyInterval(10, 500));
    set.insert(100);
    REQUIRE(!set.isEmptyInterval(10, 500));
    set.clear();
    REQUIRE(set.isEmptyInterval(10, 500));
}