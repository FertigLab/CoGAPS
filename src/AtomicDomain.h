#ifndef __COGAPS_ATOMIC_DOMAIN_H__
#define __COGAPS_ATOMIC_DOMAIN_H__

#include "data_structures/HashSets.h"
#include "math/Random.h"
#include "utils/Archive.h"

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
    ~AtomicDomain();

    // access atoms
    Atom* front();
    Atom* randomAtom(GapsRng *rng, const SmallPairedHashSetU64 &moves);
    AtomNeighborhood randomAtomWithNeighbors(GapsRng *rng);
    AtomNeighborhood randomAtomWithRightNeighbor(GapsRng *rng, const SmallPairedHashSetU64 &moves);

    Atom* getLeftNeighbor(uint64_t pos);
    Atom* getRightNeighbor(uint64_t pos);

    uint64_t randomFreePosition(GapsRng *rng,
        const std::vector<uint64_t> &possibleDeaths) const;
    uint64_t size() const;

    // these need to happen concurrently without invalidating pointers
    void erase(uint64_t pos);
    //void move(uint64_t src, uint64_t dest);

    // iterators
    std::vector<Atom*>::iterator begin() { return mAtoms.begin(); }
    std::vector<Atom*>::iterator end() { return mAtoms.end(); }
    
    bool isSorted();

    // serialization
    friend Archive& operator<<(Archive &ar, AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    // only the proposal queue can insert
    friend class ProposalQueue;
    Atom* insert(uint64_t pos, float mass);

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;

    // domain storage, sorted vector
    std::vector<Atom*> mAtoms;
};

#endif // __COGAPS_ATOMIC_DOMAIN_H__