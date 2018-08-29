#ifndef __COGAPS_ATOMIC_DOMAIN_H__
#define __COGAPS_ATOMIC_DOMAIN_H__

#include "utils/Archive.h"
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
    Atom* randomAtom(GapsRng *rng);
    AtomNeighborhood randomAtomWithNeighbors(GapsRng *rng);
    AtomNeighborhood randomAtomWithRightNeighbor(GapsRng *rng);

    uint64_t randomFreePosition(GapsRng *rng) const;
    uint64_t size() const;

    void erase(uint64_t pos);
    void insert(uint64_t pos, float mass);

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
};

#endif