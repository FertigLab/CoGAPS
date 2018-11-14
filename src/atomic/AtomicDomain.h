#ifndef __COGAPS_ATOMIC_DOMAIN_H__
#define __COGAPS_ATOMIC_DOMAIN_H__

#include "../data_structures/HashSets.h"
#include "../math/Random.h"
#include "../utils/Archive.h"

#include <vector>

// comment out to use 'new' operator for allocating atoms
#define __GAPS_USE_POOLED_ALLOCATOR__

#ifdef __GAPS_USE_POOLED_ALLOCATOR__
#include "boost/pool/object_pool.hpp"
#endif

struct Atom
{
    uint64_t pos;
    float mass;

    Atom();
    Atom(uint64_t p, float m);

    void operator=(Atom other);

    friend Archive& operator<<(Archive& ar, const Atom &a);
    friend Archive& operator>>(Archive& ar, Atom &a);
};

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

    // TODO can we have internal rng since these are always called sequentially
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
    friend Archive& operator<<(Archive &ar, const AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);
   
#ifdef GAPS_DEBUG
    std::vector<Atom*>::iterator begin();
    std::vector<Atom*>::iterator end();
    bool isSorted();
#endif

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    // only the proposal queue can insert, insert not thread-safe
    friend class ProposalQueue;
    Atom* insert(uint64_t pos, float mass);

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;

    // domain storage, sorted vector of pointers to atoms created by allocator
    std::vector<Atom*> mAtoms;
#ifdef __GAPS_USE_POOLED_ALLOCATOR__
    boost::object_pool<Atom> mAtomPool;
#else
    #warning "Not using pooled allocator for atomic"
#endif
};

#endif // __COGAPS_ATOMIC_DOMAIN_H__
