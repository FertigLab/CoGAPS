#include "AtomicDomain.h"
#include "../utils/GapsAssert.h"

#include <algorithm>
#include <limits>
#include <vector>

////////////////////////////// ATOMIC DOMAIN ///////////////////////////////////

AtomicDomain::AtomicDomain(uint64_t nBins)
{
    uint64_t binLength = std::numeric_limits<uint64_t>::max() / nBins;
    mDomainLength = binLength * nBins;
}

Atom* AtomicDomain::front()
{
    GAPS_ASSERT(size() > 0);
    return mAtomMap.begin()->second;
}

Atom* AtomicDomain::randomAtom(GapsRng *rng)
{
    GAPS_ASSERT(size() > 0);
    unsigned index = rng->uniform32(0, mAtoms.size() - 1);
    return &(mAtoms[index]);
}

uint64_t AtomicDomain::randomFreePosition(GapsRng *rng) const
{
    uint64_t pos = rng->uniform64(1, mDomainLength);
    while (mAtomMap.count(pos))
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
    unsigned index = mAtoms.size();
    mAtoms.push_back(Atom(pos, mass));
    mAtoms[index].setIndex(index);
    mAtoms[index].setIterator(mAtomMap.insert(std::pair<uint64_t, Atom*>(pos, &(mAtoms[index]))).first);

    // connect with right and left neighbors
    AtomMapType::iterator itRight(mAtoms[index].iterator());
    if (++itRight != mAtomMap.end())
    {
        mAtoms[index].setRight(itRight->second);
        itRight->second->setLeft(&(mAtoms[index]));
    }
    AtomMapType::iterator itLeft(mAtoms[index].iterator());
    if (itLeft != mAtomMap.begin())
    {
        --itLeft;
        mAtoms[index].setLeft(itLeft->second);
        itLeft->second->setRight(&(mAtoms[index]));
    }
    return &(mAtoms[index]);
}

void AtomicDomain::erase(Atom *atom)
{
    mAtomMap.erase(atom->iterator());
    if (atom->hasLeft())
    {
        atom->left()->setRight(atom->right());
    }
    if (atom->hasRight())
    {
        atom->right()->setLeft(atom->left());
    }
    
    // update the neighbors and the map
    unsigned index = atom->index();
    if (index < mAtoms.size() - 1) // we are moving the last atom
    {
        Atom *left = mAtoms.back().left();
        Atom *right = mAtoms.back().right();
        mAtoms[index] = mAtoms.back();
        mAtoms[index].setIndex(index);
        mAtoms[index].iterator()->second = &(mAtoms[index]);
        if (left != NULL)
            left->setRight(&(mAtoms[index]));
        if (right != NULL)
            right->setLeft(&(mAtoms[index]));
    }
    mAtoms.pop_back();
}

void AtomicDomain::move(Atom *atom, uint64_t newPos)
{
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
