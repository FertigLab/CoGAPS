#include <testthat.h>
#include "../testthat-tweak.h"
#include "../atomic/AtomicDomain.h"
#include "../math/Random.h"
#include "../utils/GapsPrint.h"

GapsRandomState randState(123);
GapsRng AtomicRng(&randState);



TEST_CASE("AtomicDomain populate","[atomicdomain][populate]")
{

    SECTION("Construction")
    {
        AtomicDomain domain(2);
        REQUIRE(domain.size() == 0);
    }
    SECTION("Populate")
    {
        AtomicDomain domain(2);
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
        AtomicDomain domainp(2);
        domainp.insert((uint64_t)200000,0.02); 
        //atomic coord and mass
        domainp.insert((uint64_t)400000,0.01); 
        size_t sz=domainp.size();
        Atom * del = domainp.randomAtom(&AtomicRng);
        std::cout<<"Picked position "<<del->pos()<<std::endl<<std::flush;
        REQUIRE(domainp.size() == sz);
    }
}

TEST_CASE("AtomicDomain randompick fro, empty","[atomicdomain][randompickempty]") {
    SECTION("RandomAtomChooseEmpty")
    {
        AtomicDomain domain_e(2);
        Atom * dele = domain_e.randomAtom(&AtomicRng);
        REQUIRE(domain_e.size() == 0);
    }
}

TEST_CASE("AtomicDomain depopulate","[atomicdomain][depopulate]") {
    SECTION("DePopulate")
    {
        AtomicDomain domain(2);
        domain.insert((uint64_t)200000,0.02); 
        //atomic coord and mass
        domain.insert((uint64_t)100000,0.03); 
        domain.insert((uint64_t)400000,0.01); 
        size_t sz=domain.size();
        Atom * del = domain.randomAtom(&AtomicRng);
        std::cout<<"We remove at "<<del->pos()<<std::endl<<std::flush;
        REQUIRE(domain.size() == sz);
        domain.erase(del); 
        //atomic coord
        REQUIRE(domain.size() == sz-1);
        std::cout<<"domain.size()"<<domain.size()<<std::endl<<std::flush;
    }
}

TEST_CASE("AtomicDomain depopulateempty","[atomicdomain][depopulateempty]") {
    SECTION("DePopulate")
    {
        AtomicDomain domain(2);
        Atom * del = domain.randomAtom(&AtomicRng);
        domain.erase(del); 
        //atomic coord
        REQUIRE(domain.size() == 0);
    }
}


