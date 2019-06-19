#ifndef __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__
#define __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__

#include "../data_structures/HashSets.h"
#include "../math/Random.h"
#include "../utils/Archive.h"
#include "ConcurrentAtom.h"

#include <vector>

// use a pooled allocator when creating atoms
#define __GAPS_USE_POOLED_ALLOCATOR__ 0

#if __GAPS_USE_POOLED_ALLOCATOR__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "boost/pool/object_pool.hpp"
#pragma GCC diagnostic pop
#endif

// needed for friend declarations
template <class StoragePolicy>
class AsynchronousGibbsSampler;
template <class StoragePolicy>
class SingleThreadedGibbsSampler;
class ProposalQueue;

class ConcurrentAtomicDomain
{
public:

    ConcurrentAtomicDomain(uint64_t nBins);

    // access atoms
    ConcurrentAtom* front();
    ConcurrentAtom* randomAtom(GapsRng *rng);
    ConcurrentAtomNeighborhood randomAtomWithNeighbors(GapsRng *rng);

    uint64_t randomFreePosition(GapsRng *rng) const;
    uint64_t size() const;

    // serialization
    friend Archive& operator<<(Archive &ar, const ConcurrentAtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, ConcurrentAtomicDomain &domain);
   
private:

    template <class StoragePolicy>
    friend class AsynchronousGibbsSampler;
    template <class StoragePolicy>
    friend class SingleThreadedGibbsSampler;
    friend class ProposalQueue;

    // both insert and erase will keep pointers to atoms valid, and both insert
    // and erase can be called concurrently using OpenMP threads. Neither call
    // will preserve iterators however. These protections can incur a performance
    // overhead and should only be used when needed. Friend classes are declared
    // here as a way to enforce this contract between classes.
    ConcurrentAtom* insert(uint64_t pos, float mass);
    void erase(ConcurrentAtom *atom);
    void move(ConcurrentAtom *atom, uint64_t newPos);

    ConcurrentAtomMapType mAtomMap; // sorted, used when inserting atoms to find neighbors
    std::vector<ConcurrentAtom*> mAtoms; // unsorted, used for reads

#if __GAPS_USE_POOLED_ALLOCATOR__
    boost::object_pool<ConcurrentAtom> mAtomPool;
#endif

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;
};

#endif // __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__
