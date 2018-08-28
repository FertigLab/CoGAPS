#include "catch.h"
#include "../AtomicDomain.h"
#include "../utils/GapsPrint.h"

TEST_CASE("AtomicDomain")
{
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
            REQUIRE_NOTHROW(domain.randomFreePosition());
        }
    }

    SECTION("Random Atom")
    {
        AtomicDomain domain(10);

        for (unsigned i = 0; i < 1000; ++i)
        {
            domain.insert(i, static_cast<float>(i));
            
            // single random atom
            REQUIRE(domain.randomAtom()->pos < i + 1);

            // random atom for exchange
            AtomNeighborhood hood = domain.randomAtomWithRightNeighbor();
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
            hood = domain.randomAtomWithNeighbors();
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

#if 0

// used to create aligned buckets for testing
AtomBucket* getAlignedBucket()
{
    return new AtomBucket();
}

TEST_CASE("AtomBucket")
{
    GapsRng::setSeed(12345);

    SECTION("Properties of Default Constructed Buckets")
    {
        AtomBucket *bucket = getAlignedBucket();
        REQUIRE(bucket->isEmpty());
        REQUIRE(bucket->size() == 0);
        REQUIRE(!bucket->contains(1));

        // verify internal state
        #ifdef GAPS_INTERNAL_TESTS
        REQUIRE(bucket->mOverflow == NULL);
        #endif
    }

    SECTION("Insert")
    {
        AtomBucket *bucket = getAlignedBucket();

        // insert first atom
        bucket->insert(1, 1.f);

        REQUIRE(!bucket->isEmpty());
        REQUIRE(bucket->size() == 1);

        REQUIRE(bucket->contains(1));
        REQUIRE(!bucket->contains(2));
        REQUIRE(!bucket->contains(3));

        REQUIRE(bucket->front()->pos == 1);
        REQUIRE(bucket->front()->mass == 1.f);
        REQUIRE((*bucket)[0]->pos == 1);
        REQUIRE((*bucket)[0]->mass == 1.f);

        // insert second atom
        bucket->insert(3, 3.f);

        REQUIRE(!bucket->isEmpty());
        REQUIRE(bucket->size() == 2);

        REQUIRE(bucket->contains(1));
        REQUIRE(!bucket->contains(2));
        REQUIRE(bucket->contains(3));

        REQUIRE(bucket->front()->pos == 1);
        REQUIRE(bucket->front()->mass == 1.f);
        REQUIRE((*bucket)[0]->pos == 1);
        REQUIRE((*bucket)[0]->mass == 1.f);
        REQUIRE((*bucket)[1]->pos == 3);
        REQUIRE((*bucket)[1]->mass == 3.f);

        // insert third atom
        bucket->insert(2, 2.f);

        REQUIRE(!bucket->isEmpty());
        REQUIRE(bucket->size() == 3);

        REQUIRE(bucket->contains(1));
        REQUIRE(bucket->contains(2));
        REQUIRE(bucket->contains(3));

        REQUIRE(bucket->front()->pos == 1);
        REQUIRE(bucket->front()->mass == 1.f);
        REQUIRE((*bucket)[0]->pos == 1);
        REQUIRE((*bucket)[0]->mass == 1.f);
        REQUIRE((*bucket)[1]->pos == 2);
        REQUIRE((*bucket)[1]->mass == 2.f);
        REQUIRE((*bucket)[2]->pos == 3);
        REQUIRE((*bucket)[2]->mass == 3.f);

        // verify internal state
        #ifdef GAPS_INTERNAL_TESTS
        REQUIRE(bucket->mOverflow == NULL);
        REQUIRE(bucket->back()->pos == 3);
        REQUIRE(bucket->back()->mass == 3.f);
        REQUIRE(bucket->getLeft(1)->pos == 1);
        REQUIRE(bucket->getLeft(1)->mass == 1.f);
        REQUIRE(bucket->getRight(1)->pos == 3);
        REQUIRE(bucket->getRight(1)->mass == 3.f);
        #endif
    }

    SECTION("Erase")
    {
        AtomBucket *bucket = getAlignedBucket();

        // erase one atom
        bucket->insert(1, 1.f);
        bucket->erase(1);
        REQUIRE(bucket->isEmpty());
        REQUIRE(bucket->size() == 0);
        REQUIRE(!bucket->contains(1));

        // erase two atoms
        bucket->insert(1, 1.f);
        bucket->insert(3, 3.f);

        bucket->erase(1);
        REQUIRE(!bucket->isEmpty());
        REQUIRE(bucket->size() == 1);
        REQUIRE(!bucket->contains(1));
        REQUIRE(bucket->contains(3));
        REQUIRE(bucket->front()->pos == 3);
        REQUIRE((*bucket)[0]->pos == 3);
        bucket->insert(1, 1.f);

        bucket->erase(3);
        REQUIRE(!bucket->isEmpty());
        REQUIRE(bucket->size() == 1);
        REQUIRE(bucket->contains(1));
        REQUIRE(!bucket->contains(3));
        REQUIRE(bucket->front()->pos == 1);
        REQUIRE((*bucket)[0]->pos == 1);
        bucket->insert(3, 3.f);

        bucket->erase(1);
        bucket->erase(3);
        REQUIRE(bucket->isEmpty());
        REQUIRE(bucket->size() == 0);
        REQUIRE(!bucket->contains(1));
        REQUIRE(!bucket->contains(3));
    
        // erase three atoms
        bucket->insert(1, 1.f);
        bucket->insert(3, 3.f);
        bucket->insert(2, 2.f);

        bucket->erase(2);
        REQUIRE(!bucket->isEmpty());
        REQUIRE(bucket->size() == 2);
        REQUIRE(!bucket->contains(2));
        REQUIRE(bucket->contains(1));
        REQUIRE(bucket->contains(3));
        REQUIRE(bucket->front()->pos == 1);
        REQUIRE((*bucket)[0]->pos == 1);
        REQUIRE((*bucket)[1]->pos == 3);
        bucket->insert(2, 2.f);

        bucket->erase(1);
        bucket->erase(2);
        REQUIRE(!bucket->isEmpty());
        REQUIRE(bucket->size() == 1);
        REQUIRE(!bucket->contains(1));
        REQUIRE(!bucket->contains(2));
        REQUIRE(bucket->contains(3));
        REQUIRE(bucket->front()->pos == 3);
        REQUIRE((*bucket)[0]->pos == 3);
        bucket->insert(1, 1.f);
        bucket->insert(2, 2.f);

        bucket->erase(1);
        bucket->erase(2);
        bucket->erase(3);
        REQUIRE(bucket->isEmpty());
        REQUIRE(bucket->size() == 0);
        REQUIRE(!bucket->contains(1));
        REQUIRE(!bucket->contains(2));
        REQUIRE(!bucket->contains(3));

        // verify internal state
        #ifdef GAPS_INTERNAL_TESTS
        REQUIRE(bucket->mOverflow == NULL);
        #endif
    }

    SECTION("Handling Overflow")
    {
        AtomBucket *bucket = getAlignedBucket();

        // insert 1-100 at indices 0-99 (cant use 0 for position)
        for (unsigned i = 100; i > 0; --i)
        {
            bucket->insert(i, static_cast<float>(i));
        }

        REQUIRE(!bucket->isEmpty());
        REQUIRE(bucket->size() == 100);
        REQUIRE(bucket->contains(100));
        REQUIRE(!bucket->contains(101));
        REQUIRE(bucket->front()->pos == 1);

        #ifdef GAPS_INTERNAL_TESTS
        for (unsigned i = 1; i < 99; ++i)
        {
            REQUIRE((*bucket)[i]->pos == i + 1);
            REQUIRE(bucket->contains(i + 1));
            REQUIRE(bucket->getLeft(i)->pos == i);
            REQUIRE(bucket->getRight(i)->pos == i + 2);
        }

        REQUIRE(bucket->mOverflow != NULL);
        REQUIRE(boost::alignment::is_aligned(64, bucket->mOverflow));
        REQUIRE(bucket->back()->pos == 100);
        #endif
    }

    SECTION("Special Overflow Case")
    {
        AtomBucket *bucket6 = getAlignedBucket();
        bucket6->insert(1, 1.f);
        bucket6->insert(2, 2.f);
        bucket6->insert(3, 3.f);
        bucket6->insert(4, 4.f);
        bucket6->insert(5, 5.f);
        bucket6->insert(6, 6.f);
        bucket6->insert(7, 7.f);

        REQUIRE(bucket6->size() == 7);
        bucket6->erase(7);
        REQUIRE(bucket6->size() == 6);
        bucket6->erase(6);
        REQUIRE(bucket6->size() == 5);
    }

    #ifdef GAPS_INTERNAL_TESTS
    SECTION("Traversing Adjacent Blocks")
    {
        // create buckets
        AtomBucket *first = getAlignedBucket();
        AtomBucket *second = getAlignedBucket();
        AtomBucket *third = getAlignedBucket();

        // set up order
        first->setRightAdjacentBucket(second);
        second->setRightAdjacentBucket(third);

        third->setLeftAdjacentBucket(second);
        second->setLeftAdjacentBucket(first);

        // insert into all buckets
        for (unsigned i = 10; i <= 19; ++i)
        {
            first->insert(i, static_cast<float>(i));
        }
        for (unsigned i = 20; i <= 29; ++i)
        {
            second->insert(i, static_cast<float>(i));
        }
        for (unsigned i = 30; i <= 39; ++i)
        {
            third->insert(i, static_cast<float>(i));
        }

        // check properties
        REQUIRE(!first->isEmpty());
        REQUIRE(!second->isEmpty());
        REQUIRE(!third->isEmpty());
        REQUIRE(first->size() == 10);
        REQUIRE(second->size() == 10);
        REQUIRE(third->size() == 10);

        // get neighbors across boundaries
        REQUIRE(first->getRight(9)->pos == 20);
        REQUIRE(second->getRight(9)->pos == 30);

        REQUIRE(second->getLeft(0)->pos == 19);
        REQUIRE(third->getLeft(0)->pos == 29);

        // delete middle, get neighbors again
        second->connectAdjacent();
        REQUIRE(first->getRight(9)->pos == 30);
        REQUIRE(third->getLeft(0)->pos == 19);
    }
    #endif

    SECTION("Querying Neighborhood")
    {
        // create buckets
        AtomBucket *first = getAlignedBucket();
        AtomBucket *second = getAlignedBucket();
        AtomBucket *third = getAlignedBucket();

        // insert into all buckets
        for (unsigned i = 10; i <= 19; ++i)
        {
            first->insert(i, static_cast<float>(i));
        }
        for (unsigned i = 20; i <= 29; ++i)
        {
            second->insert(i, static_cast<float>(i));
        }
        for (unsigned i = 30; i <= 39; ++i)
        {
            third->insert(i, static_cast<float>(i));
        }

        // set up order
        first->setRightAdjacentBucket(second);
        second->setRightAdjacentBucket(third);

        third->setLeftAdjacentBucket(second);
        second->setLeftAdjacentBucket(first);

        // get neighbors of border atoms
        REQUIRE(first->getNeighbors(9).left->pos == 18);
        REQUIRE(first->getNeighbors(9).center->pos == 19);
        REQUIRE(first->getNeighbors(9).right->pos == 20);
        REQUIRE(first->getRightNeighbor(9).center->pos == 19);
        REQUIRE(first->getRightNeighbor(9).right->pos == 20);

        REQUIRE(second->getNeighbors(9).left->pos == 28);
        REQUIRE(second->getNeighbors(9).center->pos == 29);
        REQUIRE(second->getNeighbors(9).right->pos == 30);
        REQUIRE(second->getRightNeighbor(9).center->pos == 29);
        REQUIRE(second->getRightNeighbor(9).right->pos == 30);

        REQUIRE(second->getNeighbors(0).left->pos == 19);
        REQUIRE(second->getNeighbors(0).center->pos == 20);
        REQUIRE(second->getNeighbors(0).right->pos == 21);
        REQUIRE(second->getRightNeighbor(0).center->pos == 20);
        REQUIRE(second->getRightNeighbor(0).right->pos == 21);

        REQUIRE(third->getNeighbors(0).left->pos == 29);
        REQUIRE(third->getNeighbors(0).center->pos == 30);
        REQUIRE(third->getNeighbors(0).right->pos == 31);
        REQUIRE(third->getRightNeighbor(0).center->pos == 30);
        REQUIRE(third->getRightNeighbor(0).right->pos == 31);

        // get neighbors of edge atoms (they shouldn't be there)
        AtomNeighborhood hood1 = first->getNeighbors(0);
        REQUIRE(!hood1.hasLeft());
        REQUIRE(hood1.hasRight());

        AtomNeighborhood hood3 = third->getNeighbors(9);
        REQUIRE(hood3.hasLeft());
        REQUIRE(!hood3.hasRight());

        #ifdef GAPS_INTERNAL_TESTS
        // delete middle, get neighbors again
        second->connectAdjacent();
        REQUIRE(first->getNeighbors(9).left->pos == 18);
        REQUIRE(first->getNeighbors(9).center->pos == 19);
        REQUIRE(first->getNeighbors(9).right->pos == 30);
        REQUIRE(first->getRightNeighbor(9).center->pos == 19);
        REQUIRE(first->getRightNeighbor(9).right->pos == 30);

        REQUIRE(third->getNeighbors(0).left->pos == 19);
        REQUIRE(third->getNeighbors(0).center->pos == 30);
        REQUIRE(third->getNeighbors(0).right->pos == 31);
        REQUIRE(third->getRightNeighbor(0).center->pos == 30);
        REQUIRE(third->getRightNeighbor(0).right->pos == 31);
        #endif
    }
}

TEST_CASE("AtomHashMap")
{
    SECTION("Default Construction")
    {
        AtomHashMap map(100);
        REQUIRE(map.size() == 0);
    }

    SECTION("Length Calculation and Hashing")
    {
        // nice size
        AtomHashMap map(100);
        map.setTotalLength(1000);
        
        #ifdef GAPS_INTERNAL_TESTS
        REQUIRE(map.mTotalLength == 1000);
        REQUIRE(map.mBucketLength == 10);

        REQUIRE(map.hash(1) == 0);
        REQUIRE(map.hash(9) == 0);
        REQUIRE(map.hash(10) == 1);
        REQUIRE(map.hash(19) == 1);
        REQUIRE(map.hash(20) == 2);
        REQUIRE(map.hash(1000) == 100);
        #endif
    }

    SECTION("Insert")
    {
        AtomHashMap map(100);
        map.setTotalLength(1000);

        map.insert(1, 1.f);
        REQUIRE(map.size() == 1);
        REQUIRE(map.contains(1));
        REQUIRE(map.front()->pos == 1);

        map.insert(2, 3.f);
        map.insert(3, 3.f);
        REQUIRE(map.size() == 3);
        REQUIRE(map.contains(2));
        REQUIRE(map.contains(3));
        REQUIRE(map.front()->pos == 1);
    }

    SECTION("Erase")
    {
        AtomHashMap map(100);
        map.setTotalLength(1000);

        map.insert(1, 1.f);
        map.insert(2, 3.f);
        map.insert(3, 3.f);

        map.erase(1);
        REQUIRE(map.size() == 2);
        REQUIRE(!map.contains(1));
        REQUIRE(map.contains(2));
        REQUIRE(map.contains(3));
        REQUIRE(map.front()->pos == 2);
    }

    #ifdef GAPS_INTERNAL_TESTS
    SECTION("Keeping Track of Full Buckets")
    {
        AtomHashMap map(100);
        map.setTotalLength(1000);

        map.insert(1, 1.f);
        REQUIRE(map.mFullBuckets.size() == 1);
        REQUIRE(map.mFullBuckets.count(0));

        map.insert(2, 2.f);
        REQUIRE(map.mFullBuckets.size() == 1);

        map.insert(10, 10.f);
        REQUIRE(map.mFullBuckets.size() == 2);
        REQUIRE(map.mFullBuckets.count(1));

        map.insert(20, 20.f);
        REQUIRE(map.mFullBuckets.size() == 3);
        REQUIRE(map.mFullBuckets.count(2));

        map.erase(10);
        REQUIRE(map.mFullBuckets.size() == 2);
        REQUIRE(map.mFullBuckets.count(0));
        REQUIRE(!map.mFullBuckets.count(1));
        REQUIRE(map.mFullBuckets.count(2));
    }
    
    SECTION("Keeping Track of Largest Bucket")
    {
        AtomHashMap map(100);
        map.setTotalLength(1000);

        for (unsigned i = 0; i < 9; ++i)
        {
            for (unsigned j = 1; j <= i + 1; ++j)
            {
                map.insert(10 * i + j, static_cast<float>(10 * i + j));
            }
            REQUIRE(map.mWhichLongestBucket == i);
            REQUIRE(map.mLongestBucketSize == i + 1);
        }

        map.erase(89);
        REQUIRE(map.mLongestBucketSize == 8);

        map.erase(88);
        map.erase(78);
        REQUIRE(map.mLongestBucketSize == 7);

        map.erase(87);
        map.erase(77);
        map.erase(67);
        REQUIRE(map.mLongestBucketSize == 6);

        map.erase(86);
        map.erase(76);
        map.erase(66);
        map.erase(56);
        REQUIRE(map.mLongestBucketSize == 5);

        map.erase(85);
        map.erase(75);
        map.erase(65);
        map.erase(55);
        map.erase(45);
        REQUIRE(map.mLongestBucketSize == 4);
    }

    SECTION("Random Index Selection")
    {
        // create has map
        AtomHashMap map(100);
        map.setTotalLength(1000);

        for (unsigned i = 0; i < 9; ++i)
        {
            for (unsigned j = 1; j <= i + 1; ++j)
            {
                map.insert(10 * i + j, static_cast<float>(10 * i + j));
            }
            REQUIRE(map.mWhichLongestBucket == i);
            REQUIRE(map.mLongestBucketSize == i + 1);
        }

        // all buckets should be selected uniformly
        unsigned counts[9] = {0};
        for (unsigned i = 0; i < 45000; ++i)
        {
            HashMapIndex index = map.getRandomIndex();
            counts[index.bucket]++;
        }

        for (unsigned i = 0; i < 9; ++i)
        {
            gaps_printf("bucket histogram: %d %d\n", i, counts[i]);
            REQUIRE(counts[i] > 950 * (i + 1));
            REQUIRE(counts[i] < 1050 * (i + 1));
        }

        // all positions should be selected uniformly
        unsigned largestBucket[9] = {0};
        for (unsigned i = 0; i < 45000; ++i)
        {
            HashMapIndex index = map.getRandomIndex();
            if (index.bucket == 8)
            {
                largestBucket[index.index]++;
            }
        }

        for (unsigned i = 0; i < 9; ++i)
        {
            gaps_printf("position histogram: %d %d\n", i, largestBucket[i]);
            REQUIRE(largestBucket[i] > 950);
            REQUIRE(largestBucket[i] < 1050);
        }
    }

    SECTION("Maintaining Adjacent Buckets")
    {
        AtomHashMap map(100);
        map.setTotalLength(1000);

        map.insert(1, 1.f);
        map.insert(9, 9.f);
        map.insert(10, 10.f);
        map.insert(19, 19.f);
        map.insert(20, 20.f);
        map.insert(29, 29.f);

        REQUIRE(map.mHashMap[0].getRight(1)->pos == 10);
        REQUIRE(map.mHashMap[1].getLeft(0)->pos == 9);

        REQUIRE(map.mHashMap[1].getRight(1)->pos == 20);
        REQUIRE(map.mHashMap[2].getLeft(0)->pos == 19);

        map.erase(10);
        map.erase(19);

        REQUIRE(map.mHashMap[0].getRight(1)->pos == 20);
        REQUIRE(map.mHashMap[2].getLeft(0)->pos == 9);
    }
    #endif

    SECTION("Random Atom Selection")
    {
        // create hash map
        AtomHashMap map(100);
        map.setTotalLength(1000);

        // average entry equal to 57
        for (unsigned i = 0; i < 9; ++i)
        {
            for (unsigned j = 1; j <= i + 1; ++j)
            {
                map.insert(10 * i + j, static_cast<float>(10 * i + j));
            }
        }

        // random selection of single atom
        float sum = 0.f;
        for (unsigned i = 0; i < 1000; ++i)
        {
            REQUIRE(map.randomAtom()->pos >= 1ull);
            REQUIRE(map.randomAtom()->pos <= 89ull);
            REQUIRE(map.randomAtom()->mass >= 1.f);
            REQUIRE(map.randomAtom()->mass <= 89.f);
            sum += map.randomAtom()->mass;
        }
        REQUIRE(sum > 0.95f * 57000.f);
        REQUIRE(sum < 1.05f * 57000.f);

        // random selection of atom and right neighbor
        sum = 0.f;
        for (unsigned i = 0; i < 1000; ++i)
        {
            AtomNeighborhood hood = map.randomAtomWithRightNeighbor();
            REQUIRE(hood.center->pos >= 1ull);
            REQUIRE(hood.center->pos <= 89ull);
            REQUIRE(hood.center->mass >= 1.f);
            REQUIRE(hood.center->mass <= 89.f);
            sum += hood.center->mass;

            REQUIRE(!hood.hasLeft());
            if (!hood.hasRight())
            {
                REQUIRE(hood.center->pos == 89);
            }
            else
            {
                REQUIRE(hood.right->pos > hood.center->pos);
                REQUIRE(hood.right->mass > hood.center->mass);
                REQUIRE(hood.right->pos - hood.center->pos <= 10);
            }
        }
        REQUIRE(sum > 0.95f * 57000.f);
        REQUIRE(sum < 1.05f * 57000.f);

        // random selection of atom and both neighbors
        sum = 0.f;
        for (unsigned i = 0; i < 1000; ++i)
        {
            AtomNeighborhood hood = map.randomAtomWithNeighbors();
            REQUIRE(hood.center->pos >= 1ull);
            REQUIRE(hood.center->pos <= 89ull);
            REQUIRE(hood.center->mass >= 1.f);
            REQUIRE(hood.center->mass <= 89.f);
            sum += hood.center->mass;

            if (!hood.hasLeft())
            {
                REQUIRE(hood.center->pos == 1);
            }
            else
            {
                REQUIRE(hood.center->pos > hood.left->pos);
                REQUIRE(hood.center->mass > hood.left->mass);
                REQUIRE(hood.center->pos - hood.left->pos <= 10);
            }

            if (!hood.hasRight())
            {
                REQUIRE(hood.center->pos == 89);
            }
            else
            {
                REQUIRE(hood.right->pos > hood.center->pos);
                REQUIRE(hood.right->mass > hood.center->mass);
                REQUIRE(hood.right->pos - hood.center->pos <= 10);
            }
        }
        REQUIRE(sum > 0.95f * 57000.f);
        REQUIRE(sum < 1.05f * 57000.f);
    }
}

#endif