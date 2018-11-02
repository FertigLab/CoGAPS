#include "catch.h"
#include "../atomic/AtomicDomain.h"
#include "../utils/GapsPrint.h"

TEST_CASE("AtomicDomain")
{
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    SECTION("Construction")
    {
        AtomicDomain domain(10);
    
        REQUIRE(domain.size() == 0);
    }

    #ifdef GAPS_INTERNAL_TESTS
    SECTION("Insert")
    {
        AtomicDomain domain(10);
        
        for (unsigned i = 0; i < 1000; ++i)
        {
            REQUIRE_NOTHROW(domain.insert(i, static_cast<float>(i)));
            REQUIRE(domain.size() == i + 1);
        }
    }

    SECTION("Erase")
    {
        AtomicDomain domain(10);
        
        for (unsigned i = 0; i < 1000; ++i)
        {
            domain.insert(i, static_cast<float>(i));
        }

        unsigned counter = domain.size();
        for (unsigned j = 0; j < 1000; j += 10)
        {
            REQUIRE_NOTHROW(domain.erase(j));
            REQUIRE(domain.size() == --counter);
        }
    }
    
    SECTION("Random Free Position")
    {
        AtomicDomain domain(10);

        for (unsigned i = 0; i < 1000; ++i)
        {
            domain.insert(i, static_cast<float>(i));
            REQUIRE_NOTHROW(domain.randomFreePosition(&rng));
        }
    }

    SECTION("Random Atom")
    {
        AtomicDomain domain(10);

        for (unsigned i = 0; i < 1000; ++i)
        {
            domain.insert(i, static_cast<float>(i));
            
            // single random atom
            REQUIRE(domain.randomAtom(&rng)->pos < i + 1);

            // random atom for exchange
            AtomNeighborhood hood = domain.randomAtomWithRightNeighbor(&rng);
            REQUIRE(hood.center->pos < i + 1);

            REQUIRE(!hood.hasLeft());
            if (hood.center->pos == i)
            {
                REQUIRE(!hood.hasRight());
            }
            else
            {
                REQUIRE(hood.hasRight());
                REQUIRE(hood.right->pos == hood.center->pos + 1);
            }

            // random atom for move
            hood = domain.randomAtomWithNeighbors(&rng);
            REQUIRE(hood.center->pos < i + 1);

            if (hood.center->pos == 0)
            {
                REQUIRE(!hood.hasLeft());
            }
            else
            {
                REQUIRE(hood.left->pos == hood.center->pos - 1);
            }

            if (hood.center->pos == i)
            {
                REQUIRE(!hood.hasRight());
            }
            else
            {
                REQUIRE(hood.hasRight());
                REQUIRE(hood.right->pos == hood.center->pos + 1);
            }
        }
    }
    #endif
}