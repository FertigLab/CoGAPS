#include "ConcurrentAtomicDomain.h"
#include "../utils/GapsAssert.h"

#include <algorithm>
#include <limits>
#include <vector>

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
    return (*mAtomMap.begin()).second;
}

Atom* ConcurrentAtomicDomain::randomAtom(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);
    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    return mAtoms[index];
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

Atom* ConcurrentAtomicDomain::insert(uint64_t pos, float mass)
{
    Atom *atom;
    #pragma omp critical(AtomicInsertOrErase)
    {
    #if __GAPS_USE_POOLED_ALLOCATOR__
        atom = mAtomPool.construct(pos, mass);
    #else
        atom = new Atom(pos, mass);
    #endif

        // insert atom into vector and map, record the iterator and index in each
        atom->setIterator(mAtomMap.insert(std::pair<uint64_t, Atom*>(pos, atom)).first);
        atom->setIndex(mAtoms.size());
        mAtoms.push_back(atom);

        // connect with right and left neighbors
        AtomMapType::iterator itRight(atom->iterator());
        if (++itRight != mAtomMap.end())
        {
            atom->setRight((*itRight).second);
            (*itRight).second->setLeft(atom);
        }
        AtomMapType::iterator itLeft(atom->iterator());
        if (itLeft != mAtomMap.begin())
        {
            --itLeft;
            atom->setLeft((*itLeft).second);
            (*itLeft).second->setRight(atom);
        }
    }
    return atom;
}

void ConcurrentAtomicDomain::erase(Atom *atom)
{
    #pragma omp critical(AtomicInsertOrErase)
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
    #if __GAPS_USE_POOLED_ALLOCATOR__
        mAtomPool.destroy(atom);
    #else
        delete atom;
    #endif
    }
}

void ConcurrentAtomicDomain::move(Atom *atom, uint64_t newPos)
{
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
    Atom temp(0, 0.f);
    uint64_t size = 0;
    ar >> domain.mDomainLength >> size;
    for (unsigned i = 0; i < size; ++i)
    {
        ar >> temp;
        domain.insert(temp.pos(), temp.mass());
    }
    return ar;
}
