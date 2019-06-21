#ifndef __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__
#define __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__

#include "../data_structures/HashSets.h"
#include "../math/Random.h"
#include "../utils/Archive.h"
#include "ConcurrentAtom.h"

#include <vector>

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
    const ConcurrentAtom* front() const;
    ConcurrentAtom* randomAtom(GapsRng *rng);
    ConcurrentAtomNeighborhood randomAtomWithNeighbors(GapsRng *rng);

    uint64_t randomFreePosition(GapsRng *rng) const;
    uint64_t size() const;

    void cacheErase(ConcurrentAtom *atom); // OpenMP thread safe
    void move(ConcurrentAtom *atom, uint64_t newPos); // OpenMP thread safe
    void flushEraseCache();

    // serialization
    friend Archive& operator<<(Archive &ar, const ConcurrentAtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, ConcurrentAtomicDomain &domain);

#ifdef GAPS_DEBUG
    bool isSorted() const;
#endif
   
private:

    // only these classes can call the non-thread safe insert and erase functions
    template <class StoragePolicy>
    friend class SingleThreadedGibbsSampler;
    friend class ProposalQueue;

    // these functions are not thread safe
    ConcurrentAtom* insert(uint64_t pos, float mass);
    void erase(ConcurrentAtom *atom);

    ConcurrentAtomMapType mAtomMap; // sorted, used when inserting atoms to find neighbors
    std::vector<ConcurrentAtom*> mAtoms; // unsorted, used for reads
    std::vector<ConcurrentAtom*> mEraseCache;

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;
};

#endif // __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__
