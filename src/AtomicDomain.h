#ifndef __COGAPS_ATOMIC_DOMAIN_H__
#define __COGAPS_ATOMIC_DOMAIN_H__

#include "Archive.h"
#include "math/Random.h"

#include <vector>

struct Atom
{
    uint64_t pos;
    float mass;

    Atom();
    Atom(uint64_t p, float m);

    void operator=(Atom other);

    friend Archive& operator<<(Archive& ar, Atom &a);
    friend Archive& operator>>(Archive& ar, Atom &a);
};

struct AtomNeighborhood
{
    Atom* center;
    Atom* left;
    Atom* right;

    AtomNeighborhood(Atom *l, Atom *c, Atom *r);

    bool hasLeft();
    bool hasRight();
};

class AtomicDomain
{
public:

    AtomicDomain(uint64_t nBins);

    // access atoms
    Atom* front();
    Atom* randomAtom();
    AtomNeighborhood randomAtomWithNeighbors();
    AtomNeighborhood randomAtomWithRightNeighbor();

    uint64_t randomFreePosition() const;
    uint64_t size() const;

    // modify domain
    void cacheInsert(uint64_t pos, float mass) const;
    void cacheErase(uint64_t pos) const;
    void flushCache();
    void resetCache(unsigned n);

    // serialization
    friend Archive& operator<<(Archive &ar, AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;

    // domain storage, sorted vector
    std::vector<Atom> mAtoms;

    // holds cache of operations
    mutable std::vector<Atom> mInsertCache;
    mutable std::vector<uint64_t> mEraseCache;

    // current index in the operation cache
    mutable unsigned mInsertCacheIndex;
    mutable unsigned mEraseCacheIndex;

    // random number generator
    mutable GapsRng mRng;

    void erase(uint64_t pos);
    void insert(uint64_t pos, float mass);

    // serialization
    friend Archive& operator<<(Archive &ar, AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);
};

#endif