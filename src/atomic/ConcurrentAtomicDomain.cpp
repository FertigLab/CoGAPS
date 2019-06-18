#include "ConcurrentAtomicDomain.h"
#include "../utils/GapsAssert.h"

#include <algorithm>
#include <limits>
#include <vector>

////////////////////////////////// HELPER //////////////////////////////////////

// used with std::lower_bound
static bool compareAtomLower(const Atom* lhs, uint64_t pos)
{
    return lhs->pos < pos;
}

// used with std::binary_search
static bool compareAtom(const Atom *lhs, const Atom *rhs)
{
    return lhs->pos < rhs->pos;
}

// check if a position in contained in a vector of atoms
static bool vecContains(const std::vector<Atom*> &vec, uint64_t pos)
{
    Atom temp(pos, 0.f);
    return std::binary_search(vec.begin(), vec.end(), &temp, compareAtom);
}

////////////////////////////////// Atom ////////////////////////////////////////

Atom::Atom() : pos(0), mass(0.f) {}

Atom::Atom(uint64_t p, float m) : pos(p), mass(m) {}

void Atom::operator=(Atom other)
{
    pos = other.pos;
    mass = other.mass;
}

Archive& operator<<(Archive &ar, const Atom &a)
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

AtomNeighborhood::AtomNeighborhood()
    : center(NULL), left(NULL), right(NULL)
{}

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

#if __GAPS_USE_POOLED_ALLOCATOR__
ConcurrentAtomicDomain::ConcurrentAtomicDomain(uint64_t nBins)
    : mAtomPool(16384, 0),
#else
ConcurrentAtomicDomain::ConcurrentAtomicDomain(uint64_t nBins)
#endif
{
    uint64_t binLength = std::numeric_limits<uint64_t>::max() / nBins;
    mDomainLength = binLength * nBins;
}

Atom* ConcurrentAtomicDomain::front()
{
    GAPS_ASSERT(size() > 0);

    return mAtoms.front();
}

Atom* ConcurrentAtomicDomain::randomAtom(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);

    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    return mAtoms[index];
}

AtomNeighborhood ConcurrentAtomicDomain::randomAtomWithNeighbors(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);

    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    Atom* left = (index == 0) ? NULL : mAtoms[index - 1];
    Atom* right = (index == mAtoms.size() - 1) ? NULL : mAtoms[index + 1];
    return AtomNeighborhood(left, mAtoms[index], right);
}

AtomNeighborhood ConcurrentAtomicDomain::randomAtomWithRightNeighbor(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);

    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    Atom* right = (index == mAtoms.size() - 1) ? NULL : mAtoms[index + 1];
    return AtomNeighborhood(NULL, mAtoms[index], right);
}

uint64_t ConcurrentAtomicDomain::randomFreePosition(GapsRng *rng) const
{
    uint64_t pos = rng->uniform64(1, mDomainLength);
    while (vecContains(mAtoms, pos))
    {
        pos = rng->uniform64(1, mDomainLength);
    }
    return pos;
}

uint64_t ConcurrentAtomicDomain::size() const
{
    return mAtoms.size();
}

Atom* ConcurrentAtomicDomain::insert(uint64_t pos, float mass)
{
    Atom *insertedAtom = NULL;
    #pragma omp critical(AtomicInsertOrErase)
    {
        std::vector<Atom*>::iterator it;
        GAPS_ASSERT(!vecContains(mAtoms, pos));
        it = std::lower_bound(mAtoms.begin(), mAtoms.end(), pos, compareAtomLower);
    #if __GAPS_USE_POOLED_ALLOCATOR__
        it = mAtoms.insert(it, mAtomPool.construct(pos, mass));
    #else
        it = mAtoms.insert(it, new Atom(pos, mass));
    #endif
        GAPS_ASSERT(vecContains(mAtoms, pos));
        insertedAtom = *it;
    }
    return insertedAtom;
}

void ConcurrentAtomicDomain::erase(uint64_t pos)
{
    #pragma omp critical(AtomicInsertOrErase)
    {
        GAPS_ASSERT(size() > 0);
        GAPS_ASSERT_MSG(vecContains(mAtoms, pos), pos);
        std::vector<Atom*>::iterator it = std::lower_bound(mAtoms.begin(),
            mAtoms.end(), pos, compareAtomLower);
    #if __GAPS_USE_POOLED_ALLOCATOR__
        mAtomPool.destroy(*it);
    #else
        delete *it;
    #endif
        mAtoms.erase(it);
        GAPS_ASSERT(!vecContains(mAtoms, pos));
    }
}

Archive& operator<<(Archive &ar, const ConcurrentAtomicDomain &domain)
{
    ar << domain.mDomainLength << domain.mAtoms.size();
    for (unsigned i = 0; i < domain.mAtoms.size(); ++i)
    {
        ar << *(domain.mAtoms[i]);
    }
    return ar;
}

Archive& operator>>(Archive &ar, ConcurrentAtomicDomain &domain)
{
    Atom temp;
    uint64_t size = 0;
    ar >> domain.mDomainLength >> size;
    for (unsigned i = 0; i < size; ++i)
    {
        ar >> temp;
        domain.insert(temp.pos, temp.mass);
    }
    return ar;
}

////////////////////// ATOMIC DOMAIN DEBUG FUNCTIONS ///////////////////////////

#ifdef GAPS_DEBUG

std::vector<Atom*>::iterator ConcurrentAtomicDomain::begin()
{
    return mAtoms.begin();
}

std::vector<Atom*>::iterator ConcurrentAtomicDomain::end()
{
    return mAtoms.end();
}

// used in debug mode to check if vector is always sorted
bool ConcurrentAtomicDomain::isSorted()
{
    for (unsigned i = 1; i < mAtoms.size(); ++i)
    {
        if (mAtoms[i]->pos <= mAtoms[i-1]->pos)
        {
            gaps_printf("unsorted\n%llu\n%llu\n", mAtoms[i-1]->pos, mAtoms[i]->pos);
            return false;
        }
    }
    return true;
}

#endif

