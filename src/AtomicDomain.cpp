#include "AtomicDomain.h"
#include "Random.h"

#include <stdint.h>
#include <utility>

// O(1)
Atom AtomicDomain::front() const
{
    return mAtoms[mAtomPositions.begin()->second];
}

// O(1)
Atom AtomicDomain::randomAtom() const
{
    uint64_t num = gaps::random::uniform64(0, mAtoms.size() - 1);
    return mAtoms[num];
}

// O(logN) - keep hash of positions to fix, need this O(1)
uint64_t AtomicDomain::randomFreePosition() const
{
    uint64_t pos = 0;
    do
    {
        pos = gaps::random::uniform64();
    } while (mAtomPositions.count(pos) > 0); // count is O(logN)
    return pos;
}

// O(logN)
void AtomicDomain::insert(uint64_t pos, float mass)
{
    // insert position into map
    std::map<uint64_t, uint64_t>::iterator iter, iterLeft, iterRight;
    iter = mAtomPositions.insert(std::pair<uint64_t, uint64_t>(pos, mAtoms.size())).first;
    iterLeft = iter;
    iterRight = iter;

    // find neighbors
    Atom atom(pos, mass);
    if (iter != mAtomPositions.begin())
    {
        --iterLeft;
        atom.left = &(mAtoms[iterLeft->second]);
    }
    if (iter != mAtomPositions.end())
    {
        ++iterRight;
        atom.right = &(mAtoms[iterRight->second]);
    } 

    // add atom to vector
    mAtoms.push_back(atom);
}

// O(logN)
void AtomicDomain::erase(uint64_t pos)
{
    // get vector index of this atom and erase it
    uint64_t index = mAtomPositions.at(pos);
    mAtomPositions.erase(pos);

    // update key of object about to be moved (last one doesn't need to move)
    if (index < mAtoms.size() - 1)
    {
        mAtomPositions.erase(mAtoms.back().pos);
        mAtomPositions.insert(std::pair<uint64_t, uint64_t>(mAtoms.back().pos,
            index));
    }

    // update neighbors
    if (mAtoms[index].left)
    {
        mAtoms[index].left->right = mAtoms[index].right;
    }
    if (mAtoms[index].right)
    {
        mAtoms[index].right->left = mAtoms[index].left;
    }

    // delete atom from vector in O(1)
    mAtoms[index] = mAtoms.back();
    mAtoms.pop_back();
}

// O(logN)
void AtomicDomain::updateMass(uint64_t pos, float newMass)
{
    uint64_t index = mAtomPositions.at(pos);
    mAtoms[index].mass = newMass;
}
