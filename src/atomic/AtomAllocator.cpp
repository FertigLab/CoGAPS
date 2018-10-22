#include "AtomAllocator.h"

////////////////////////////////// Atom ////////////////////////////////////////

Atom::Atom() : pos(0), mass(0.f) {}

Atom::Atom(uint64_t p, float m) : pos(p), mass(m) {}

void Atom::operator=(Atom other)
{
    pos = other.pos;
    mass = other.mass;
}

Archive& operator<<(Archive &ar, Atom &a)
{
    ar << a.pos << a.mass;
    return ar;
}

Archive& operator>>(Archive &ar, Atom &a)
{
    ar >> a.pos >> a.mass;
    return ar;
}

////////////////////////////// IndexFlagSet ////////////////////////////////////

// all bits start as free
IndexFlagSet::IndexFlagSet()
{
    for (unsigned i = 0; i < NUM_INDEX_CHUNKS; ++i)
    {
        parts[i] = -1;
    }
}

// find first set bit
// TODO this is compiler dependent
unsigned IndexFlagSet::getFirstFree() const
{
    GAPS_ASSERT(isAnyFree());

    for (unsigned i = 0; i < NUM_INDEX_CHUNKS; ++i)
    {
        if (parts[i] != 0)
        {
            return __builtin_ffsll(parts[i]) + 64 * i - 1;
        }
    }
    return 0; // ERROR IF REACHED
}

bool IndexFlagSet::isAnyFree() const
{
#if (NUM_INDEX_CHUNKS == 8)
    return (parts[0] | parts[1] | parts[2] | parts[3] | parts[4] | parts[5]
        | parts[6] | parts[7]) != 0;
#else
    #warning "non-standard size for atom pool"
    for (unsigned i = 0; i < NUM_INDEX_CHUNKS; ++i)
    {
        if (parts[i] != 0)
        {
            return true;
        }
    }
    return false;
#endif
}

// set position n to 0
void IndexFlagSet::set(uint64_t n)
{
    unsigned i = n / 64;
    parts[i] ^= (1ull << (n - 64 * i));
}

// set position n to 1
void IndexFlagSet::release(uint64_t n)
{
    unsigned i = n / 64;
    parts[i] |= (1ull << (n - 64 * i));
}

//////////////////////////////// AtomPool //////////////////////////////////////

AtomPool::AtomPool()
{
    mPool = new Atom[POOL_SIZE];
}

AtomPool::~AtomPool()
{
    delete[] mPool;
}

Atom* AtomPool::alloc()
{
    unsigned n = mIndexFlags.getFirstFree();
    mIndexFlags.set(n);
    Atom *a = &(mPool[n]);
    a->poolIndex = n;
    return a;
}

void AtomPool::free(Atom* a)
{
    mIndexFlags.release(a->poolIndex);
}

bool AtomPool::depleted() const
{
    return !mIndexFlags.isAnyFree();
}

////////////////////////////// AtomAllocator ///////////////////////////////////

// for debugging
#define __USE_CUSTOM_ALLOCATOR__ 1

AtomAllocator::AtomAllocator()
{
#if __USE_CUSTOM_ALLOCATOR__
    mIndex = 0;
    mPools.push_back(new AtomPool());
#endif
}

AtomAllocator::~AtomAllocator()
{
#if __USE_CUSTOM_ALLOCATOR__
    std::vector<AtomPool*>::iterator it = mPools.begin();
    for (; it != mPools.end(); ++it)
    {
        delete *it;
    }
#endif
}

Atom* AtomAllocator::createAtom(uint64_t pos, float mass)
{
#if __USE_CUSTOM_ALLOCATOR__
    GAPS_ASSERT(mPools.size() < 65536);

    if (mPools[mIndex]->depleted() && mIndex == mPools.size() - 1)
    {
        mPools.push_back(new AtomPool());
        mIndex = 0; // loop back around all pools before getting to new one
    }

    while (mPools[mIndex]->depleted())
    {
        ++mIndex;
    }

    Atom *a = mPools[mIndex]->alloc();
    a->allocatorIndex = mIndex;
    a->pos = pos;
    a->mass = mass;
    return a;
#else
    return new Atom(pos, mass);
#endif
}

void AtomAllocator::destroyAtom(Atom *a)
{
#if __USE_CUSTOM_ALLOCATOR__
    mPools[a->allocatorIndex]->free(a);
#else
    delete a;
#endif
}
