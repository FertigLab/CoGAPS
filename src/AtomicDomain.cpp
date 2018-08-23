#include "AtomicDomain.h"
#include "GapsAssert.h"

#include <vector>

////////////////////////////////// ATOM ////////////////////////////////////////

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

//////////////////////////// ATOM NEIGHBORHOOD /////////////////////////////////


AtomNeighborhood::AtomNeighborhood(Atom *l, Atom *c, Atom *r)
    : center(c), left(l), right(r)
{}

bool AtomNeighborhood::hasLeft()
{
    return left != NULL;
}

bool AtomNeighborhood::hasRight()
{
    return right != NULL;
}

////////////////////////////// ATOMIC DOMAIN ///////////////////////////////////

AtomicDomain::AtomicDomain(unsigned nBins)
{
    uint64_t binLength = std::numeric_limits<uint64_t>::max() / nBins;
    mDomainLength = binLength * nBins;
}

Atom* AtomicDomain::front()
{
    GAPS_ASSERT(size() > 0);
    return &(mAtoms.front());
}

Atom* AtomicDomain::randomAtom()
{
    GAPS_ASSERT(size() > 0);
    unsigned index = mRng.uniform32(0, mAtoms.size() - 1);
    return &(mAtoms[index]);
}

AtomNeighborhood AtomicDomain::randomAtomWithNeighbors()
{
    GAPS_ASSERT(size() > 0);
    unsigned index = mRng.uniform32(0, mAtoms.size() - 1);
    Atom* left = (index == 0) ? NULL : &(mAtoms[index - 1]);
    Atom* right = (index == mAtoms.size() - 1) ? NULL : &(mAtoms[index + 1]);
    return AtomNeighborhood(left, &(mAtoms[index]), right);
}

AtomNeighborhood AtomicDomain::randomAtomWithRightNeighbor()
{
    GAPS_ASSERT(size() > 0);
    unsigned index = mRng.uniform32(0, mAtoms.size() - 1);
    Atom* right = (index == mAtoms.size() - 1) ? NULL : &(mAtoms[index + 1]);
    return AtomNeighborhood(NULL, &(mAtoms[index]), right);
}

static bool compareAtomLower(const Atom &lhs, uint64_t pos)
{
    return lhs.pos < pos;
}

static bool compareAtomUpper(uint64_t pos, const Atom &lhs)
{
    return lhs.pos < pos;
}

static bool compareAtom(const Atom &lhs, const Atom &rhs)
{
    return lhs.pos < rhs.pos;
}

static bool vecContains(const std::vector<Atom> &vec, uint64_t pos)
{
    Atom temp(pos, 0.f);
    return std::binary_search(vec.begin(), vec.end(), temp, compareAtom);
}

uint64_t AtomicDomain::randomFreePosition() const
{
    uint64_t pos = mRng.uniform64(0, mDomainLength);
    while (vecContains(mAtoms, pos))
    {
        pos = mRng.uniform64(0, mDomainLength);
    } 
    return pos;
}

uint64_t AtomicDomain::size() const
{
    return mAtoms.size();
}

void AtomicDomain::cacheInsert(uint64_t pos, float mass) const
{
    unsigned ndx = 0;
    #pragma omp critical(atomicInsert)
    {
        ndx = mInsertCacheIndex++;
    }
    mInsertCache[ndx] = Atom(pos, mass);
}

void AtomicDomain::cacheErase(uint64_t pos) const
{
    unsigned ndx = 0;
    #pragma omp critical(atomicErase)
    {
        ndx = mEraseCacheIndex++;
    }
    mEraseCache[ndx] = pos;
}

void AtomicDomain::resetCache(unsigned n)
{
    mInsertCacheIndex = 0;
    mEraseCacheIndex = 0;
    mInsertCache.resize(n);
    mEraseCache.resize(n);
}

void AtomicDomain::erase(uint64_t pos)
{
    GAPS_ASSERT(size() > 0);
    std::vector<Atom>::iterator it;
    it = std::lower_bound(mAtoms.begin(), mAtoms.end(), pos, compareAtomLower);
    mAtoms.erase(it);
}

void AtomicDomain::insert(uint64_t pos, float mass)
{
    std::vector<Atom>::iterator it;
    it = std::lower_bound(mAtoms.begin(), mAtoms.end(), pos, compareAtomLower);
    mAtoms.insert(it, Atom(pos, mass));
}

void AtomicDomain::flushCache()
{
    for (unsigned i = 0; i < mEraseCacheIndex; ++i)
    {
        erase(mEraseCache[i]);
    }

    for (unsigned i = 0; i < mInsertCacheIndex; ++i)
    {
        insert(mInsertCache[i].pos, mInsertCache[i].mass);
    }

    mInsertCache.clear();
    mEraseCache.clear();
}

Archive& operator<<(Archive &ar, AtomicDomain &domain)
{
    ar << domain.mDomainLength << domain.mRng << domain.mAtoms.size();
    
    for (unsigned i = 0; i < domain.mAtoms.size(); ++i)
    {
        ar << domain.mAtoms[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, AtomicDomain &domain)
{
    Atom temp;
    unsigned size = 0;
    ar >> domain.mDomainLength >> domain.mRng >> size;
    for (unsigned i = 0; i < size; ++i)
    {
        ar >> temp;
        domain.insert(temp.pos, temp.mass);
    }
    return ar;
}
