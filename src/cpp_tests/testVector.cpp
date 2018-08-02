#include "catch.h"
#include "../data_structures/Vector.h"

TEST_CASE("Test Vector.h")
{
    Vector v1(100);
    Vector v2(std::vector<float>(100.f, 4));

    SECTION("Test Construction")
    {
        REQUIRE(v1.size() == 100);
        REQUIRE(v1[0] == 0.f);

        REQUIRE(v2.size() == 4);
        REQUIRE(v2[0] == 100.f);
    }
    
    SECTION("Test Concatenation")
    {
        v1.concat(v2);
        REUQIRE(v1.size() == 104);
    }
}