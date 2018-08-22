#ifndef __COGAPS_ATOMIC_DOMAIN_H__
#define __COGAPS_ATOMIC_DOMAIN_H__

#include "Archive.h"
#include "math/Random.h"

#include <vector>
#include <set>

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

class AtomBucket
{
public:

    AtomBucket();

    unsigned size() const;
    bool isEmpty() const;
    bool contains(uint64_t pos) const;

    Atom* front();
    Atom* operator[](unsigned index);

    void insert(uint64_t pos, float mass);
    void erase(uint64_t pos);
        
    AtomNeighborhood getNeighbors(unsigned index);
    AtomNeighborhood getRightNeighbor(unsigned index);

    void setRightAdjacentBucket(AtomBucket *bucket);
    void setLeftAdjacentBucket(AtomBucket *bucket);

    AtomBucket& operator=(const AtomBucket& other);

    unsigned getIndex(uint64_t pos);
    Atom* getLeft(unsigned index);
    Atom* getRight(unsigned index);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    Atom mBucket;
    bool mFull;
    
    AtomBucket *mOverflow;
    AtomBucket *mPrev;
    AtomBucket *mNext;

    void eraseFront();
    void connectAdjacent();

    Atom* back();
    
    friend Archive& operator<<(Archive& ar, AtomBucket &b);
    friend Archive& operator>>(Archive& ar, AtomBucket &b);
};

struct HashMapIndex
{
    unsigned bucket;
    unsigned index;

    HashMapIndex(unsigned b, unsigned i) : bucket(b), index(i) {}
};

// hash map of chained AtomBuckets
// data structure that holds atoms
// accessing happens much more than insert/erase so more time is spent
// up front (sorting) to save time on access
// note that atoms can never occupy position 0
class AtomHashMap
{
public:

    // required that length is a multiple of expectedNumAtoms
    AtomHashMap(unsigned expectedNumAtoms);
    void setTotalLength(uint64_t length);

    // TODO while moving the atoms cannot change the atom order, it can
    // move an atom from one bucket to another - this moved atom must be
    // the front atom moving to the back of the prev bucket or the back
    // atom moving to the front of the next atom - effectively creating
    // a reorder event - it may even need allocation

    Atom* front();
    Atom* randomAtom();
    AtomNeighborhood randomAtomWithNeighbors();
    AtomNeighborhood randomAtomWithRightNeighbor();

    bool contains(uint64_t pos) const;
    unsigned size() const;

    void insert(uint64_t pos, float mass);
    void erase(uint64_t pos);
    void move(uint64_t src, uint64_t dest, float mass); // cannot reorder

    Atom* getLeft(uint64_t pos);
    Atom* getRight(uint64_t pos);

    bool hasLeft(uint64_t pos);
    bool hasRight(uint64_t pos);

    void updateMass(uint64_t pos, float newMass);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    unsigned mSize;
    unsigned mWhichLongestBucket;
    unsigned mLongestBucketSize;

    unsigned mExpectedNumAtoms;

    uint64_t mBucketLength;
    uint64_t mTotalLength;

    std::vector<AtomBucket> mHashMap;
    std::set<unsigned> mFullBuckets; // TODO use IntFixedHashSet

    // random number generator
    mutable GapsRng mRng;

    unsigned hash(uint64_t pos) const;
    HashMapIndex getRandomIndex() const;

    friend Archive& operator<<(Archive& ar, AtomHashMap &h);
    friend Archive& operator>>(Archive& ar, AtomHashMap &h);
};

struct MoveOperation
{
    uint64_t src;
    uint64_t dest;
    float mass;

    MoveOperation() : src(0), dest(0), mass(0.f) {}

    MoveOperation(uint64_t s, uint64_t d, float m) :
        src(s), dest(d), mass(m)
    {}
};

class AtomicDomain
{
public:

    AtomicDomain(unsigned nBins);

    void setDomainSize(uint64_t size);

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
    void cacheMove(uint64_t src, uint64_t dest, float mass) const;
    void flushCache();
    void resetCache(unsigned n);

    // serialization
    friend Archive& operator<<(Archive &ar, AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);

private:

    // size of atomic domain to ensure all bins are equal length
    uint64_t mDomainLength;

    // domain storage, specialized hash map
    mutable AtomHashMap mAtoms;

    // holds cache of operations
    mutable std::vector<Atom> mInsertCache;
    mutable std::vector<uint64_t> mEraseCache;
    mutable std::vector<MoveOperation> mMoveCache;

    // current index in the operation cache
    mutable unsigned mInsertCacheIndex;
    mutable unsigned mEraseCacheIndex;
    mutable unsigned mMoveCacheIndex;

    // random number generator
    mutable GapsRng mRng;

    // serialization
    friend Archive& operator<<(Archive &ar, AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);
};

#endif