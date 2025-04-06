#include "AtomicDomain.h"
#include "../utils/GapsAssert.h"
#include "../utils/Archive.h"
#include "../math/Random.h"

#include <algorithm>
#include <limits>
#include <vector>

////////////////////////////// ATOMIC DOMAIN ///////////////////////////////////

AtomicDomain::AtomicDomain(uint64_t nBins)
{
    //nBins is the naumber of matrix elements to be fitted (nrows(A)*ncols(A)+nrows(P)*ncols(P)) in our run;
    //each element has a band in atomic space; the max possible atomic coord is the max uit_64 and
    //we want the bands to be equal, so see below 
    uint64_t binLength = std::numeric_limits<uint64_t>::max() / nBins;
    mDomainLength = binLength * nBins;
}

Atom* AtomicDomain::front()
{
    GAPS_ASSERT(size() > 0);
    return &(mAtoms[mAtomMap.begin()->second]);
    //we take the first element in atom map, take the value (we use key->value notation in comments),
    //and index the atom storage by the value and return the reference to the atom
}

Atom* AtomicDomain::randomAtom(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);
    size_t index = rng->uniform64(0, mAtoms.size() - 1);
    return &(mAtoms[index]);
}

AtomNeighborhood AtomicDomain::randomAtomWithNeighbors(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);
    size_t index = rng->uniform64(0, mAtoms.size() - 1);
    Atom *center = &(mAtoms[index]);
    Atom *left = center->hasLeft() ? &(mAtoms[center->leftIndex()]) : NULL;
    Atom *right = center->hasRight() ? &(mAtoms[center->rightIndex()]) : NULL;
    return AtomNeighborhood(left, center, right);
}

uint64_t AtomicDomain::randomFreePosition(GapsRng *rng) const
{
    uint64_t pos = rng->uniform64(1, mDomainLength);
    while (mAtomMap.count(pos) != 0u)
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
    size_t index = mAtoms.size();
    //it will be the index of the atom if inserted
    auto insertresult=mAtomMap.insert(std::pair<uint64_t, size_t>(pos, index));
    if (!insertresult.second) {
        //already exists - so the second ib the returned pair is false
        //and the first is the map iterator to the existing pair
        index=insertresult.first->second;
        return &(mAtoms[index]);
    }
    //if we are here, we did de novo insertion
    mAtoms.push_back(Atom(pos, mass));
    mAtoms[index].setIndex(index);
    mAtoms[index].setIterator(insertresult.first);

    // connect with right and left neighbors
    AtomMapType::iterator itRight(mAtoms[index].iterator());
    if (++itRight != mAtomMap.end())
    {
        mAtoms[index].setRightIndex(itRight->second);
        mAtoms[itRight->second].setLeftIndex(index);
    }
    AtomMapType::iterator itLeft(mAtoms[index].iterator());
    if (itLeft != mAtomMap.begin())
    {
        --itLeft;
        mAtoms[index].setLeftIndex(itLeft->second);
        mAtoms[itLeft->second].setRightIndex(index);
    }
    return &(mAtoms[index]);
}

void AtomicDomain::erase(Atom *atom)
{
    mAtomMap.erase(atom->iterator());
    if (atom->hasLeft())
    {
        mAtoms[atom->leftIndex()].setRightIndex(atom->rightIndex());
    }
    if (atom->hasRight())
    {
        mAtoms[atom->rightIndex()].setLeftIndex(atom->leftIndex());
    }
    
    // update the neighbors and the map
    size_t index = atom->index();
    if (index < mAtoms.size() - 1) // we are moving the last atom
    {
        int leftIndex = mAtoms.back().leftIndex();
        int rightIndex = mAtoms.back().rightIndex();
        mAtoms[index] = mAtoms.back();
        mAtoms[index].setIndex(index);
        mAtoms[index].iterator()->second = index;
        if (leftIndex >= 0)
        {
            mAtoms[leftIndex].setRightIndex(index);
        }
        if (rightIndex >= 0)
        {
            mAtoms[rightIndex].setLeftIndex(index);
        }
    }
    mAtoms.pop_back();
}

void AtomicDomain::move(Atom *atom, uint64_t newPos)
{
    GAPS_ASSERT(newPos > (atom->hasLeft() ? mAtoms[atom->leftIndex()].pos() : 0));
    GAPS_ASSERT(newPos < (atom->hasRight() ? mAtoms[atom->rightIndex()].pos() : mDomainLength));
    atom->updatePos(newPos);
    mAtomMap.updateKey(atom->iterator(), newPos);
}

Archive& operator<<(Archive &ar, const AtomicDomain &domain)
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
