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

    uint64_t randomFreePosition(GapsRng *rng) const;
    uint64_t size() const;

    // serialization
    friend Archive& operator<<(Archive &ar, const AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);
   
private:

    template <class StoragePolicy>
    friend class SingleThreadedGibbsSampler;

    // both insert and erase will invalidate pointers to atoms and neither
    // of these functions should be considered thread-safe in any way. This
    // class can provide increased performance when thread-safety and 
    // concurrency is not needed. Friend classes are declared here as a way
    // to enforce this contract between classes.
    Atom* insert(uint64_t pos, float mass);
    void erase(Atom *atom);
    void move(Atom *atom, uint64_t newPos);

    AtomMapType mAtomMap; // sorted, used when inserting atoms to find neighbors
    std::vector<Atom> mAtoms; // unsorted, used for reads

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;
};

#endif // __COGAPS_ATOMIC_DOMAIN_H__
