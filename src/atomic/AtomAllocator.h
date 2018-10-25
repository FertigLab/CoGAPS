#ifndef __COGAPS_ATOM_ALLOCATOR_H__
#define __COGAPS_ATOM_ALLOCATOR_H__

#include "../utils/Archive.h"

#include <stdint.h>
#include <vector>

// NOTE: can only have 65536 pools with a 16bit index, so POOL_SIZE * 65536
// is the maximum number of atoms CoGAPS can have - don't change this unless
// you know what you are doing
#define POOL_SIZE 512 // allows for 33.5 million atoms
#define NUM_INDEX_CHUNKS POOL_SIZE / 64

class AtomAllocator;
class AtomPool;

struct Atom
{
public:

    uint64_t pos;
    float mass;

    Atom();
    Atom(uint64_t p, float m);

    void operator=(Atom other);

    friend Archive& operator<<(Archive& ar, const Atom &a);
    friend Archive& operator>>(Archive& ar, Atom &a);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    // used by the allocator
    friend class AtomAllocator;
    friend class AtomPool;
    uint16_t poolIndex;
    uint16_t allocatorIndex;
};

// stores flags for POOL_SIZE indices
struct IndexFlagSet
{
public:

    IndexFlagSet();

    unsigned getFirstFree() const;
    bool isAnyFree() const;
    void set(uint64_t n);
    void release(uint64_t n);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    uint64_t parts[NUM_INDEX_CHUNKS];
};

class AtomPool
{
public:
    
    AtomPool();
    ~AtomPool();

    Atom* alloc();
    void free(Atom* a);
    bool depleted() const;

#ifndef GAPS_INTERNAL_TESTS
private:
#endif
    
    Atom* mPool;
    IndexFlagSet mIndexFlags;

    AtomPool(const AtomPool &pool); // = delete (no c++11)
    AtomPool& operator=(const AtomPool &pool); // = delete (no c++11)
};

class AtomAllocator
{
public:

    AtomAllocator();
    ~AtomAllocator();

    Atom* createAtom(uint64_t pos, float mass);
    void destroyAtom(Atom *a);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    std::vector<AtomPool*> mPools;
    unsigned mIndex;

    AtomAllocator(const AtomAllocator &pool); // = delete (no c++11)
    AtomAllocator& operator=(const AtomAllocator &pool); // = delete (no c++11)
};

#endif