#include "AtomicDomain.h"
#include "utils/GapsAssert.h"

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

AtomicDomain::~AtomicDomain()
{
    for (unsigned i = 0; i < mAtoms.size(); ++i)
    {
        delete mAtoms[i];
    }
}

Atom* AtomicDomain::front()
{
    GAPS_ASSERT(size() > 0);

    return mAtoms.front();
}

Atom* AtomicDomain::randomAtom(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);
    GAPS_ASSERT(isSorted());

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

uint64_t AtomicDomain::randomFreePosition(GapsRng *rng,
const std::vector<uint64_t> &possibleDeaths) const
{
    uint64_t pos = rng->uniform64(1, mDomainLength);
    while (vecContains(mAtoms, pos))
    {
        if (vecContains(possibleDeaths, pos))
        {
            return 0; // might actually be a free position
        }
        pos = rng->uniform64(1, mDomainLength);
    }
    return pos;
}

uint64_t AtomicDomain::size() const
{
    return mAtoms.size();
}

void AtomicDomain::erase(uint64_t pos)
{
    Atom *a = NULL;
    #pragma omp critical(AtomicInsertOrErase)
    {
        GAPS_ASSERT(size() > 0);
        GAPS_ASSERT_MSG(vecContains(mAtoms, pos), pos);

        std::vector<Atom*>::iterator it = std::lower_bound(mAtoms.begin(),
            mAtoms.end(), pos, compareAtomLower);
        Atom *a = *it;
        mAtoms.erase(it);

        GAPS_ASSERT(!vecContains(mAtoms, pos));
    }
    delete a;
}

void AtomicDomain::move(uint64_t src, uint64_t dest)
{
    #pragma omp critical(AtomicInsertOrErase)
    {
        GAPS_ASSERT(size() > 0);
        GAPS_ASSERT_MSG(vecContains(mAtoms, src), src);
        GAPS_ASSERT(!vecContains(mAtoms, dest));

        std::vector<Atom*>::iterator it = std::lower_bound(mAtoms.begin(),
            mAtoms.end(), src, compareAtomLower);
        unsigned ndx = std::distance(mAtoms.begin(), it);
        while (ndx + 1 < mAtoms.size() && dest > mAtoms[ndx + 1]->pos)
        {
            Atom* temp = mAtoms[ndx];
            mAtoms[ndx] = mAtoms[ndx + 1];
            mAtoms[ndx + 1] = temp;
            ++ndx;
        }
        while (ndx > 0 && dest < mAtoms[ndx - 1]->pos)
        {
            Atom* temp = mAtoms[ndx];
            mAtoms[ndx] = mAtoms[ndx - 1];
            mAtoms[ndx - 1] = temp;
            --ndx;
        }
        mAtoms[ndx]->pos = dest;

        GAPS_ASSERT(!vecContains(mAtoms, src));
        GAPS_ASSERT(vecContains(mAtoms, dest));
    }
}

Atom* AtomicDomain::insert(uint64_t pos, float mass)
{
    Atom *newAtom = new Atom(pos, mass);
    std::vector<Atom*>::iterator it;
    #pragma omp critical(AtomicInsertOrErase)
    {
        GAPS_ASSERT(!vecContains(mAtoms, pos));

        it = std::lower_bound(mAtoms.begin(), mAtoms.end(), pos, compareAtomLower);
        it = mAtoms.insert(it, newAtom);

        GAPS_ASSERT(vecContains(mAtoms, pos));
    }
    return *it;
}

Archive& operator<<(Archive &ar, AtomicDomain &domain)
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
