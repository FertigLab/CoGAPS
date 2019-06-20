#include "ConcurrentAtomicDomain.h"
#include "../utils/GapsAssert.h"

#include <algorithm>
#include <limits>

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

ConcurrentAtom* ConcurrentAtomicDomain::insert(uint64_t pos, float mass)
{
    ConcurrentAtom *atom;
    #pragma omp critical(AtomicInsertOrErase)
    {
        // insert atom into vector and map, record the iterator and index in each
        atom = new ConcurrentAtom(pos, mass);
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
    }
    return atom;
}

void ConcurrentAtomicDomain::erase(ConcurrentAtom *atom)
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
        delete atom;
    }
}

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
