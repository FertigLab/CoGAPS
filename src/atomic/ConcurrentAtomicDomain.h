#ifndef __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__
#define __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__

#include "ConcurrentAtom.h"

#include <vector>

template <class StoragePolicy>
class SingleThreadedGibbsSampler;

struct ConcurrentAtom;
class Archive;
class GapsRng;
class ProposalQueue;

class ConcurrentAtomicDomain
{
public:
    explicit ConcurrentAtomicDomain(uint64_t nBins);
    ConcurrentAtom* front();
    const ConcurrentAtom* front() const;
    ConcurrentAtom* randomAtom(GapsRng *rng);
    ConcurrentAtomNeighborhood randomAtomWithNeighbors(GapsRng *rng);
    uint64_t randomFreePosition(GapsRng *rng) const;
    uint64_t size() const;
    void cacheErase(ConcurrentAtom *atom); // OpenMP thread safe
    void move(ConcurrentAtom *atom, uint64_t newPos); // OpenMP thread safe
    void flushEraseCache();
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
    std::vector<ConcurrentAtom*> mAtoms; // unsorted, used for random selection of atoms
    std::vector<ConcurrentAtom*> mEraseCache;
    uint64_t mDomainLength; // size of atomic domain to ensure all bins are equal length
};

#endif // __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__
