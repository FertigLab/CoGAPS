#ifndef __COGAPS_ATOMIC_DOMAIN_H__
#define __COGAPS_ATOMIC_DOMAIN_H__

#include "AtomAllocator.h"
#include "data_structures/HashSets.h"
#include "math/Random.h"
#include "math/SIMD.h"
#include "utils/Archive.h"

#include <vector>

struct AtomNeighborhood
{
    Atom* center;
    Atom* left;
    Atom* right;

    AtomNeighborhood();
    AtomNeighborhood(Atom *l, Atom *c, Atom *r);

    bool hasLeft();
    bool hasRight();
};

class ProposalQueue; // needed for friend

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

    // this needs to happen concurrently without invalidating pointers
    void erase(uint64_t pos);

    // serialization
    friend Archive& operator<<(Archive &ar, AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);
   
#ifdef GAPS_DEBUG
    std::vector<Atom*>::iterator begin();
    std::vector<Atom*>::iterator end();
    bool isSorted();
#endif

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    // only the proposal queue can insert
    friend class ProposalQueue;
    Atom* insert(uint64_t pos, float mass);

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;

    // domain storage, sorted vector of pointers to atoms created by allocator
    std::vector<Atom*> mAtoms;
    AtomAllocator mAllocator;
};

#endif // __COGAPS_ATOMIC_DOMAIN_H__
