#include "AtomicDomain.h"
#include "../utils/GapsAssert.h"

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

// check if a position in contained in a vector of positions
static bool vecContains(const std::vector<uint64_t> &vec, uint64_t pos)
{
    return std::binary_search(vec.begin(), vec.end(), pos);
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

AtomicDomain::AtomicDomain(uint64_t nBins)
{
    uint64_t binLength = std::numeric_limits<uint64_t>::max() / nBins;
    mDomainLength = binLength * nBins;
}

Atom* AtomicDomain::front()
{
    GAPS_ASSERT(size() > 0);

    return mAtoms.front();
}

Atom* AtomicDomain::randomAtom(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);

    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    return mAtoms[index];
}

AtomNeighborhood AtomicDomain::randomAtomWithNeighbors(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);

    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    Atom* left = (index == 0) ? NULL : mAtoms[index - 1];
    Atom* right = (index == mAtoms.size() - 1) ? NULL : mAtoms[index + 1];
    return AtomNeighborhood(left, mAtoms[index], right);
}

AtomNeighborhood AtomicDomain::randomAtomWithRightNeighbor(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);

    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    Atom* right = (index == mAtoms.size() - 1) ? NULL : mAtoms[index + 1];
    return AtomNeighborhood(NULL, mAtoms[index], right);
}

uint64_t AtomicDomain::randomFreePosition(GapsRng *rng) const
{
    uint64_t pos = rng->uniform64(1, mDomainLength);
    while (vecContains(mAtoms, pos))
    {
        pos = rng->uniform64(1, mDomainLength);
    }
    return pos;
}

uint64_t AtomicDomain::size() const
{
    return mAtoms.size();
}

Atom* AtomicDomain::insert(uint64_t pos, float mass)
{
    GAPS_ASSERT(!vecContains(mAtoms, pos));

    std::vector<Atom*>::iterator it;
    it = std::lower_bound(mAtoms.begin(), mAtoms.end(), pos, compareAtomLower);
    it = mAtoms.insert(it, mAllocator.createAtom(pos, mass));

    GAPS_ASSERT(vecContains(mAtoms, pos));
    return *it;
}

void AtomicDomain::erase(uint64_t pos)
{
    #pragma omp critical(AtomicInsertOrErase)
    {
        GAPS_ASSERT(size() > 0);
        GAPS_ASSERT_MSG(vecContains(mAtoms, pos), pos);

        std::vector<Atom*>::iterator it = std::lower_bound(mAtoms.begin(),
            mAtoms.end(), pos, compareAtomLower);
        mAllocator.destroyAtom(*it);
        mAtoms.erase(it);

        GAPS_ASSERT(!vecContains(mAtoms, pos));
    }
}

Archive& operator<<(Archive &ar, const AtomicDomain &domain)
{
    ar << domain.mDomainLength << domain.mAtoms.size();
    
    for (unsigned i = 0; i < domain.mAtoms.size(); ++i)
    {
        ar << *(domain.mAtoms[i]);
    }
    return ar;
}

Archive& operator>>(Archive &ar, AtomicDomain &domain)
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

std::vector<Atom*>::iterator AtomicDomain::begin()
{
    return mAtoms.begin();
}

std::vector<Atom*>::iterator AtomicDomain::end()
{
    return mAtoms.end();
}

// used in debug mode to check if vector is always sorted
bool AtomicDomain::isSorted()
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

