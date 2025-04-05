#include <testthat.h>
#include "../testthat-tweak.h"
#include "../atomic/AtomicDomain.h"
#include "../math/Random.h"
#include "../utils/GapsPrint.h"

TEST_CASE("AtomicDomain-basic","[atomicdomain][basic]")
{
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    AtomicDomain domain(2);
    SECTION("Construction")
    {
        REQUIRE(domain.size() == 0);
    }
    SECTION("Populate")
    {
        domain.insert((uint64_t)100000,0.01); 
        //atomic coord and mass
        REQUIRE(domain.size() == 1);
        domain.insert((uint64_t)200000,0.02); 
        //atomic coord and mass
        REQUIRE(domain.size() == 2);
        domain.insert((uint64_t)400000,0.01); 
        //atomic coord and mass
        REQUIRE(domain.size() == 3);
        domain.insert((uint64_t)400000,0.05); 
        //atomic coord and mass
        REQUIRE(domain.size() == 3);
    }
    SECTION("DePopulate")
    {
        size_t sz=domain.size();
        //domain.erase((uint64_t)100000); 
        //atomic coord
        //REQUIRE(domain.size() == sz-1);
    }
}

