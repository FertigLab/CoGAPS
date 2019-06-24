#ifndef __COGAPS_ATOMIC_DOMAIN_H__
#define __COGAPS_ATOMIC_DOMAIN_H__

#include "Atom.h"

#include <vector>
#include <cstddef>

template <class StoragePolicy>
class SingleThreadedGibbsSampler;

class Archive;
class GapsRng;

class AtomicDomain
{
public:
    explicit AtomicDomain(uint64_t nBins);
    Atom* front();
    Atom* randomAtom(GapsRng *rng);
    AtomNeighborhood randomAtomWithNeighbors(GapsRng *rng);
    uint64_t randomFreePosition(GapsRng *rng) const;
    uint64_t size() const;
    friend Archive& operator<<(Archive &ar, const AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);
private:
    // both insert and erase will invalidate pointers to atoms and neither
    // of these functions should be considered thread-safe in any way. This
    // class can provide increased performance when thread-safety and 
    // concurrency is not needed. Friend classes are declared here as a way
    // to enforce this contract between classes.
    template <class StoragePolicy>
    friend class SingleThreadedGibbsSampler;

    Atom* insert(uint64_t pos, float mass);
    void erase(Atom *atom);
    void move(Atom *atom, uint64_t newPos);

    AtomMapType mAtomMap; // sorted, used when inserting atoms to find neighbors
    std::vector<Atom> mAtoms; // unsorted, used for reads
    uint64_t mDomainLength; // size of atomic domain to ensure all bins are equal length
};

#endif // __COGAPS_ATOMIC_DOMAIN_H__
