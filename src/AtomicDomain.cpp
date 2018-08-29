#include "AtomicDomain.h"
#include "utils/GapsAssert.h"

#include <vector>

////////////////////////////////// HELPER //////////////////////////////////////

// used with std::lower_bound
static bool compareAtomLower(const Atom &lhs, uint64_t pos)
{
    return lhs.pos < pos;
}

// used with std::binary_search
static bool compareAtom(const Atom &lhs, const Atom &rhs)
{
    return lhs.pos < rhs.pos;
}

// check if a position in contained in a vector of atoms
static bool vecContains(const std::vector<Atom> &vec, uint64_t pos)
{
    Atom temp(pos, 0.f);
    return std::binary_search(vec.begin(), vec.end(), temp, compareAtom);
}

// used in debug mode to check if vector is always sorted
static bool isSorted(const std::vector<Atom> &vec)
{
    for (unsigned i = 1; i < vec.size(); ++i)
    {
        if (vec[i].pos <= vec[i-1].pos)
        {
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

    return &(mAtoms.front());
}

Atom* AtomicDomain::randomAtom(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);
    GAPS_ASSERT(isSorted(mAtoms));

    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    return &(mAtoms[index]);
}

AtomNeighborhood AtomicDomain::randomAtomWithNeighbors(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);

    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    Atom* left = (index == 0) ? NULL : &(mAtoms[index - 1]);
    Atom* right = (index == mAtoms.size() - 1) ? NULL : &(mAtoms[index + 1]);
    return AtomNeighborhood(left, &(mAtoms[index]), right);
}

AtomNeighborhood AtomicDomain::randomAtomWithRightNeighbor(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);

    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    Atom* right = (index == mAtoms.size() - 1) ? NULL : &(mAtoms[index + 1]);
    return AtomNeighborhood(NULL, &(mAtoms[index]), right);
}

uint64_t AtomicDomain::randomFreePosition(GapsRng *rng) const
{
    uint64_t pos = rng->uniform64(0, mDomainLength);
    while (vecContains(mAtoms, pos))
    {
        pos = rng->uniform64(0, mDomainLength);
    } 
    return pos;
}

uint64_t AtomicDomain::size() const
{
    return mAtoms.size();
}

void AtomicDomain::erase(uint64_t pos)
{
    GAPS_ASSERT(size() > 0);
    GAPS_ASSERT(vecContains(mAtoms, pos));

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

Archive& operator<<(Archive &ar, AtomicDomain &domain)
{
    ar << domain.mDomainLength << domain.mAtoms.size();
    
    for (unsigned i = 0; i < domain.mAtoms.size(); ++i)
    {
        ar << domain.mAtoms[i];
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
