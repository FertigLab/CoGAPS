#include "ConcurrentAtomicDomain.h"
#include "../math/Random.h"
#include "../utils/Archive.h"
#include "../utils/GapsAssert.h"

#include <algorithm>
#include <limits>

static bool compareAtoms(ConcurrentAtom *a1, ConcurrentAtom *a2)
{
    return a1->pos() < a2->pos();
}

ConcurrentAtomicDomain::ConcurrentAtomicDomain(uint64_t nBins)
{
    uint64_t binLength = std::numeric_limits<uint64_t>::max() / nBins;
    mDomainLength = binLength * nBins;
}

ConcurrentAtom* ConcurrentAtomicDomain::front()
{
    GAPS_ASSERT(size() > 0);
    return (*mAtomMap.begin()).second;
}

const ConcurrentAtom* ConcurrentAtomicDomain::front() const
{
    GAPS_ASSERT(size() > 0);
    return (*mAtomMap.begin()).second;
}

ConcurrentAtom* ConcurrentAtomicDomain::randomAtom(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);
    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    return mAtoms[index];
}

ConcurrentAtomNeighborhood ConcurrentAtomicDomain::randomAtomWithNeighbors(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);
    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    return ConcurrentAtomNeighborhood(mAtoms[index]->left(), mAtoms[index], mAtoms[index]->right());
}

uint64_t ConcurrentAtomicDomain::randomFreePosition(GapsRng *rng) const
{
    uint64_t pos = rng->uniform64(1, mDomainLength);
    while (mAtomMap.count(pos))
    {
        pos = rng->uniform64(1, mDomainLength);
    }
    return pos;
}

uint64_t ConcurrentAtomicDomain::size() const
{
    return mAtoms.size();
}

// safe to call concurrently from OpenMP threads
void ConcurrentAtomicDomain::cacheErase(ConcurrentAtom *atom)
{
    #pragma omp critical(AtomicInsertOrErase)
    {
        mEraseCache.push_back(atom);
    }
}

// not thread safe
void ConcurrentAtomicDomain::flushEraseCache()
{
    std::sort(mEraseCache.begin(), mEraseCache.end(), compareAtoms);
    for (unsigned i = 0; i < mEraseCache.size(); ++i)
    {
        erase(mEraseCache[i]);
    }
    mEraseCache.clear();
}

// not thread safe
ConcurrentAtom* ConcurrentAtomicDomain::insert(uint64_t pos, float mass)
{
    // insert atom into vector and map, record the iterator and index in each
    ConcurrentAtom *atom = new ConcurrentAtom(pos, mass);
    atom->setIterator(mAtomMap.insert(std::pair<uint64_t, ConcurrentAtom*>(pos, atom)).first);
    atom->setIndex(mAtoms.size());
    mAtoms.push_back(atom);

    // connect with right and left neighbors
    ConcurrentAtomMapType::iterator itRight(atom->iterator());
    if (++itRight != mAtomMap.end())
    {
        atom->setRight((*itRight).second);
        (*itRight).second->setLeft(atom);
    }
    ConcurrentAtomMapType::iterator itLeft(atom->iterator());
    if (itLeft != mAtomMap.begin())
    {
        --itLeft;
        atom->setLeft((*itLeft).second);
        (*itLeft).second->setRight(atom);
    }
    return atom;
}

// not thread safe
void ConcurrentAtomicDomain::erase(ConcurrentAtom *atom)
{
    mAtomMap.erase(atom->iterator());
    mAtoms[atom->index()] = mAtoms.back();
    mAtoms[atom->index()]->setIndex(atom->index());
    mAtoms.pop_back();
    if (atom->hasLeft())
    {
        atom->left()->setRight(atom->right());
    }
    if (atom->hasRight())
    {
        atom->right()->setLeft(atom->left());
    }
    delete atom;
}

// safe to call concurrently from OpenMP threads
void ConcurrentAtomicDomain::move(ConcurrentAtom *atom, uint64_t newPos)
{
    GAPS_ASSERT(newPos > (atom->hasLeft() ? atom->left()->pos() : 0));
    GAPS_ASSERT(newPos < (atom->hasRight() ? atom->right()->pos() : mDomainLength));
    atom->updatePos(newPos);
    mAtomMap.updateKey(atom->iterator(), newPos);
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
    ConcurrentAtom temp(0, 0.f);
    uint64_t size = 0;
    ar >> domain.mDomainLength >> size;
    for (unsigned i = 0; i < size; ++i)
    {
        ar >> temp;
        domain.insert(temp.pos(), temp.mass());
    }
    return ar;
}

#ifdef GAPS_DEBUG
bool ConcurrentAtomicDomain::isSorted() const
{
    if (size() == 0)
        return true;
    const ConcurrentAtom *atom = front();
    bool valid = !atom->hasLeft();
    while (atom->hasRight())
    {
        valid = valid && (atom->right()->pos() > atom->pos());
        atom = atom->right();
        valid = valid && (atom->left()->pos() < atom->pos());
    }
    return valid;
}
#endif