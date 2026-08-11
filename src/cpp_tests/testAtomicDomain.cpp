#include <testthat.h>
#include <Rcpp.h>
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

TEST_CASE("AtomicDomain move","[atomicdomain][move]")
{
    AtomicDomain domain(100);
    domain.insert((uint64_t)200000,0.01); 
    //atomic coord and mass
    domain.insert((uint64_t)100000,0.02); 
    //atomic coord and mass
    domain.insert((uint64_t)400000,0.03); 
    //atomic coord and mass
    domain.insert((uint64_t)300000,0.04); 
    //domain.move(a2p,30000);
    REQUIRE(domain.size() == 4);
    
    SECTION("Test move")
    {
        domain.move(domain.storedAtom(3),310000);
        REQUIRE(domain.storedAtom(3)->pos() == 310000);
        auto iter=domain.front_it();
        REQUIRE(iter->first==100000);
        REQUIRE(iter->second==1);
        iter++;
        REQUIRE(iter->first==200000);
        REQUIRE(iter->second==0);
        iter++;
        REQUIRE(iter->first==310000);
        REQUIRE(iter->second==3);
        iter++;
        REQUIRE(iter->first==400000);
        REQUIRE(iter->second==2);
    }
    
}

TEST_CASE("AtomicDomain structure","[atomicdomain][structure]")
{
    AtomicDomain domain(100);
    domain.insert((uint64_t)200000,0.01); 
    //atomic coord and mass
    domain.insert((uint64_t)100000,0.02); 
    //atomic coord and mass
    domain.insert((uint64_t)400000,0.03); 
    //atomic coord and mass
    domain.insert((uint64_t)300000,0.04); 
    //domain.move(a2p,30000);
    REQUIRE(domain.size() == 4);
    
    SECTION("Test structure")
    {
        REQUIRE(domain.front()->pos() == 100000);
        REQUIRE(domain.front()->index() == 1);
        REQUIRE(domain.front()->hasLeft() == false);
        REQUIRE(domain.front()->hasRight() == true);
        REQUIRE(domain.front()->rightIndex() == 0); 
        
        REQUIRE(domain.storedAtom(0)->pos() == 200000);
        REQUIRE(domain.storedAtom(0)->index() == 0);
        REQUIRE(domain.storedAtom(0)->hasLeft() == true);
        REQUIRE(domain.storedAtom(0)->hasRight() == true);
        REQUIRE(domain.storedAtom(0)->leftIndex() == 1);
        REQUIRE(domain.storedAtom(0)->rightIndex() == 3);  
        
        REQUIRE(domain.storedAtom(1)->pos() == 100000);
        REQUIRE(domain.storedAtom(1)->index() == 1);
        REQUIRE(domain.storedAtom(1)->hasLeft() == false);
        REQUIRE(domain.storedAtom(1)->hasRight() == true);
        //REQUIRE(domain.storedAtom(1)->leftIndex() == 100);
        REQUIRE(domain.storedAtom(1)->rightIndex() == 0);
        
        REQUIRE(domain.storedAtom(2)->pos() == 400000);
        REQUIRE(domain.storedAtom(2)->index() == 2);
        REQUIRE(domain.storedAtom(2)->hasLeft() == true);
        REQUIRE(domain.storedAtom(2)->hasRight() == false);
        REQUIRE(domain.storedAtom(2)->leftIndex() == 3);
        //REQUIRE(domain.storedAtom(2)->rightIndex() == 100);
        
        REQUIRE(domain.storedAtom(3)->pos() == 300000);
        REQUIRE(domain.storedAtom(3)->index() == 3);
        REQUIRE(domain.storedAtom(3)->hasLeft() == true);
        REQUIRE(domain.storedAtom(3)->hasRight() == true);
        REQUIRE(domain.storedAtom(3)->leftIndex() == 0);
        REQUIRE(domain.storedAtom(3)->rightIndex() == 2);         
    }
    SECTION("Test map structure")
    {
        auto iter=domain.front_it();
        REQUIRE(iter->first==100000);
        REQUIRE(iter->second==1);
        iter++;
        REQUIRE(iter->first==200000);
        REQUIRE(iter->second==0);
        iter++;
        REQUIRE(iter->first==300000);
        REQUIRE(iter->second==3);
        iter++;
        REQUIRE(iter->first==400000);
        REQUIRE(iter->second==2);
    }
    SECTION("Test move")
    {
        domain.move(domain.storedAtom(3),310000);
        REQUIRE(domain.storedAtom(3)->pos() == 310000);
        auto iter=domain.front_it();
        REQUIRE(iter->first==100000);
        REQUIRE(iter->second==1);
        iter++;
        REQUIRE(iter->first==200000);
        REQUIRE(iter->second==0);
        iter++;
        REQUIRE(iter->first==310000);
        REQUIRE(iter->second==3);
        iter++;
        REQUIRE(iter->first==400000);
        REQUIRE(iter->second==2);
    }
    
}


TEST_CASE("AtomicDomain depopulate","[atomicdomain][depopulate]") {
    AtomicDomain domain(100);
    domain.insert((uint64_t)200000,0.01); 
    //atomic coord and mass
    domain.insert((uint64_t)100000,0.02); 
    //atomic coord and mass
    domain.insert((uint64_t)400000,0.03); 
    //atomic coord and mass
    domain.insert((uint64_t)300000,0.04); 
    //domain.move(a2p,30000);
    REQUIRE(domain.size() == 4);
    SECTION("DePopulate")
    {   
        domain.erase(domain.storedAtom(0));
        
        REQUIRE(domain.front()->pos() == 100000);
        REQUIRE(domain.front()->index() == 1);
        REQUIRE(domain.front()->hasLeft() == false);
        REQUIRE(domain.front()->hasRight() == true);
        REQUIRE(domain.front()->rightIndex() == 0); 
        
        REQUIRE(domain.storedAtom(0)->pos() == 300000);
        REQUIRE(domain.storedAtom(0)->index() == 0);
        REQUIRE(domain.storedAtom(0)->hasLeft() == true);
        REQUIRE(domain.storedAtom(0)->hasRight() == true);
        REQUIRE(domain.storedAtom(0)->leftIndex() == 1);
        REQUIRE(domain.storedAtom(0)->rightIndex() == 2);  
        
        REQUIRE(domain.storedAtom(1)->pos() == 100000);
        REQUIRE(domain.storedAtom(1)->index() == 1);
        REQUIRE(domain.storedAtom(1)->hasLeft() == false);
        REQUIRE(domain.storedAtom(1)->hasRight() == true);
        //REQUIRE(domain.storedAtom(1)->leftIndex() == 100);
        REQUIRE(domain.storedAtom(1)->rightIndex() == 0);
        
        REQUIRE(domain.storedAtom(2)->pos() == 400000);
        REQUIRE(domain.storedAtom(2)->index() == 2);
        REQUIRE(domain.storedAtom(2)->hasLeft() == true);
        REQUIRE(domain.storedAtom(2)->hasRight() == false);
        REQUIRE(domain.storedAtom(2)->leftIndex() == 0);
        //REQUIRE(domain.storedAtom(2)->rightIndex() == 100);
    }
    SECTION("Test map structure")
    {
        //in new section, we get a new copy of domain
        //to be sure,
        REQUIRE(domain.size()==4);
        //erase
        domain.erase(domain.storedAtom(0));
        auto iter=domain.front_it();
        REQUIRE(iter->first==100000);
        REQUIRE(iter->second==1);
        iter++;
        REQUIRE(iter->first==300000);
        REQUIRE(iter->second==0);
        iter++;
        REQUIRE(iter->first==400000);
        REQUIRE(iter->second==2);
    }
}

TEST_CASE("AtomicDomain move then erase", "[atomicdomain][movethenerase]")
{
    // Regression test: after move(), the atom's internal map iterator must be
    // updated so that a subsequent erase() on the same atom does not use a
    // stale iterator (stale iterator caused UB / map corruption).
    AtomicDomain domain(100);
    domain.insert((uint64_t)200000, 0.01); // storedAtom(0)
    domain.insert((uint64_t)100000, 0.02); // storedAtom(1)
    domain.insert((uint64_t)400000, 0.03); // storedAtom(2)
    domain.insert((uint64_t)300000, 0.04); // storedAtom(3)
    REQUIRE(domain.size() == 4);

    Atom *atom = domain.storedAtom(3); // atom at 300000
    REQUIRE(atom->pos() == 300000);
    domain.move(atom, 310000);
    REQUIRE(atom->pos() == 310000);

    // Erase the moved atom — must use the updated iterator, not the stale one
    domain.erase(atom);
    REQUIRE(domain.size() == 3);

    // 310000 must be absent; map must contain exactly 100000, 200000, 400000
    auto iter = domain.front_it();
    REQUIRE(iter->first == 100000);
    iter++;
    REQUIRE(iter->first == 200000);
    iter++;
    REQUIRE(iter->first == 400000);
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

TEST_CASE("AtomicDomain randompick from empty","[atomicdomain][randompickempty]") {
    SECTION("RandomAtomChooseEmpty")
    {
        AtomicDomain domain_e(2);
        REQUIRE(domain_e.size() == 0);
        REQUIRE_THROWS(domain_e.randomAtom(&AtomicRng));
    }
}



TEST_CASE("AtomicDomain depopulateempty","[atomicdomain][depopulateempty]") {
    SECTION("DePopulate")
    {   
        AtomicDomain domain(2);
        REQUIRE(domain.size() == 0);
        Atom * del;
        REQUIRE_THROWS(del=domain.randomAtom(&AtomicRng));
        REQUIRE_THROWS(domain.erase(del)); 
    }
}


