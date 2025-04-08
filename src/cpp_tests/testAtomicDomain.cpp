#include <testthat.h>
#include "../testthat-tweak.h"
#include "../atomic/AtomicDomain.h"
#include "../math/Random.h"
#include "../utils/GapsPrint.h"

GapsRandomState randState(123);
GapsRng AtomicRng(&randState);

AtomicDomain domain(2);


TEST_CASE("AtomicDomain populate","[atomicdomain][populate]")
{

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
}

TEST_CASE("AtomicDomain randompick","[atomicdomain][randompick]") {
    SECTION("RandomAtomChoose")
    {
        size_t sz=domain.size();
        Atom * del = domain.randomAtom(&AtomicRng);
        std::cout<<"Picked position "<<del->pos()<<std::endl<<std::flush;
        REQUIRE(domain.size() == sz);
    }
}

TEST_CASE("AtomicDomain depopulate","[atomicdomain][depopulate]") {
    SECTION("DePopulate")
    {
        size_t sz=domain.size();
        Atom * del = domain.randomAtom(&AtomicRng);
        std::cout<<"We remove at "<<del->pos()<<std::endl<<std::flush;
        REQUIRE(domain.size() == sz);
        domain.erase(del); 
        //atomic coord
        REQUIRE(domain.size() == sz-1);
        std::cout<<"domain.size()"<<domain.size()<<std::endl<<std::flush;
        //domain.erase(domain.randomAtom(&AtomicRng)); 
        //atomic coord
        //REQUIRE(domain.size() == sz-2);
        
        //domain.erase(domain.randomAtom(&AtomicRng)); 
        //atomic coord
        //REQUIRE(domain.size() == sz-3);
        
    }
}

