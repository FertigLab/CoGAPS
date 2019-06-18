#ifndef __COGAPS_ATOMIC_DOMAIN_H__
#define __COGAPS_ATOMIC_DOMAIN_H__

#include "../data_structures/HashSets.h"
#include "../math/Random.h"
#include "../utils/Archive.h"
#include "Atom.h"

#include <vector>

// needed for friend declarations
template <class StoragePolicy>
class SingleThreadedGibbsSampler;

class AtomicDomain
{
public:

    AtomicDomain(uint64_t nBins);

    // access atoms
    Atom* front();
    Atom* randomAtom(GapsRng *rng);
    AtomNeighborhood randomAtomWithNeighbors(GapsRng *rng);
    AtomNeighborhood randomAtomWithRightNeighbor(GapsRng *rng);

    uint64_t randomFreePosition(GapsRng *rng) const;
    uint64_t size() const;

    // serialization
    friend Archive& operator<<(Archive &ar, const AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);
   
#ifdef GAPS_DEBUG
    std::vector<Atom*>::iterator begin();
    std::vector<Atom*>::iterator end();
    bool isSorted();
#endif

private:

    // both insert and erase will invalidate pointers to atoms and neither
    // of these functions should be considered thread-safe in any way. This
    // class can provide increased performance when thread-safety and 
    // concurrency is not needed. Friend classes are declared here as a way
    // to enforce this contract between classes.
    template <class StoragePolicy>
    friend class SingleThreadedGibbsSampler;
    Atom* insert(uint64_t pos, float mass);
    void erase(uint64_t pos);

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;

    // domain storage, sorted vector of pointers to atoms created by allocator
    std::vector<Atom> mAtoms;
};

#endif // __COGAPS_ATOMIC_DOMAIN_H__
