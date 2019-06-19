#ifndef __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__
#define __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__

#include "../data_structures/HashSets.h"
#include "../math/Random.h"
#include "../utils/Archive.h"
#include "Atom.h"

#include <vector>

// use a pooled allocator when creating atoms
#define __GAPS_USE_POOLED_ALLOCATOR__ 0

#if __GAPS_USE_POOLED_ALLOCATOR__
#include "boost/pool/object_pool.hpp"
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
    Atom* front();
    Atom* randomAtom(GapsRng *rng);

    uint64_t randomFreePosition(GapsRng *rng) const;
    uint64_t size() const;

    // serialization
    friend Archive& operator<<(Archive &ar, const ConcurrentAtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, ConcurrentAtomicDomain &domain);
   
private:

    // both insert and erase will keep pointers to atoms valid, and both insert
    // and erase can be called concurrently using OpenMP threads. Neither call
    // will preserve iterators however. These protections can incur a performance
    // overhead and should only be used when needed. Friend classes are declared
    // here as a way to enforce this contract between classes.
    template <class StoragePolicy>
    friend class AsynchronousGibbsSampler;
    template <class StoragePolicy>
    friend class SingleThreadedGibbsSampler;
    friend class ProposalQueue;
    Atom* insert(uint64_t pos, float mass);
    void erase(Atom *atom);
    void move(Atom *atom, uint64_t newPos);

    AtomMapType mAtomMap; // sorted, used when inserting atoms to find neighbors
    std::vector<Atom*> mAtoms; // unsorted, used for reads

#if __GAPS_USE_POOLED_ALLOCATOR__
    boost::object_pool<Atom> mAtomPool;
#endif

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;
};

#endif // __COGAPS_CONCURRENT_ATOMIC_DOMAIN_H__
