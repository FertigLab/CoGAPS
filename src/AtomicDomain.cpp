#include "GapsAssert.h"
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
    GAPS_ASSERT(mAtoms.size() > 0);
    uint64_t num = gaps::random::uniform64(0, mAtoms.size() - 1);
    GAPS_ASSERT(num >= 0);
    GAPS_ASSERT(num < mAtoms.size());
    return mAtoms[num];
}

// Average Case O(1)
uint64_t AtomicDomain::randomFreePosition() const
{
    uint64_t pos = 0;
    do
    {
        pos = gaps::random::uniform64();
    } while (mUsedPositions.count(pos) > 0); // hash map => count is O(l)
    return pos;
}

// O(1)
uint64_t AtomicDomain::size() const
{
    return mAtoms.size();
}

// O(1)
Atom& AtomicDomain::left(const Atom &atom)
{
    GAPS_ASSERT(hasLeft(atom));
    return mAtoms[atom.leftNdx - 1];
}

// O(1)
Atom& AtomicDomain::right(const Atom &atom)
{
    GAPS_ASSERT(hasRight(atom));
    return mAtoms[atom.rightNdx - 1];
}

// O(1)
const Atom& AtomicDomain::left(const Atom &atom) const
{
    GAPS_ASSERT(hasLeft(atom));
    return mAtoms[atom.leftNdx - 1];
}

// O(1)
const Atom& AtomicDomain::right(const Atom &atom) const
{
    GAPS_ASSERT(hasRight(atom));
    return mAtoms[atom.rightNdx - 1];
}

// O(1)
bool AtomicDomain::hasLeft(const Atom &atom) const
{
    return atom.leftNdx != 0;
}

// O(1)
bool AtomicDomain::hasRight(const Atom &atom) const
{
    return atom.rightNdx != 0;
}

// O(logN)
void AtomicDomain::insert(uint64_t pos, float mass)
{
    // insert position into map
    GAPS_ASSERT(mAtomPositions.find(pos) == mAtomPositions.end());
    std::map<uint64_t, uint64_t>::iterator iter, iterLeft, iterRight;
    iter = mAtomPositions.insert(std::pair<uint64_t, uint64_t>(pos, size())).first;
    iterLeft = iter;
    iterRight = iter;

    // find neighbors
    Atom atom(pos, mass);
    if (iter != mAtomPositions.begin())
    {
        --iterLeft;
        atom.leftNdx = iterLeft->second + 1;
        left(atom).rightNdx = size() + 1;
    }
    if (++iter != mAtomPositions.end())
    {
        ++iterRight;
        atom.rightNdx = iterRight->second + 1;
        right(atom).leftNdx = size() + 1;
    } 

    // add atom to vector
    GAPS_ASSERT(atom.rightNdx <= size());
    GAPS_ASSERT(atom.leftNdx <= size());
    mAtoms.push_back(atom);
    mUsedPositions.insert(pos);
    GAPS_ASSERT(mAtoms.size() == mAtomPositions.size());
    GAPS_ASSERT(mAtoms.size() == mUsedPositions.size());
    GAPS_ASSERT(atom.rightNdx != size() + 1);
    GAPS_ASSERT(atom.leftNdx != size() + 1);
}

// O(logN)
// erase directly from position map and used positions hash set
// swap with last atom in vector, pop off the back
void AtomicDomain::erase(uint64_t pos)
{
    // get vector index of this atom and erase it
    GAPS_ASSERT(mAtomPositions.find(pos) != mAtomPositions.end());
    uint64_t index = mAtomPositions.at(pos);

    // connect neighbors of atom to be deleted
    if (hasLeft(mAtoms[index]))
    {
        left(mAtoms[index]).rightNdx = mAtoms[index].rightNdx;
    }
    if (hasRight(mAtoms[index]))
    {
        right(mAtoms[index]).leftNdx = mAtoms[index].leftNdx;
    }

    // replace with atom from back
    if (index < size() - 1)
    {
        mAtoms[index] = mAtoms.back();

        // update position in map
        mAtomPositions.erase(mAtoms[index].pos);
        mAtomPositions.insert(std::pair<uint64_t, uint64_t>(mAtoms[index].pos,
            index));
    
        // update moved atom's neighbors
        if (hasLeft(mAtoms[index]))
        {
            left(mAtoms[index]).rightNdx = index + 1;
        }
        if (hasRight(mAtoms[index]))
        {
            right(mAtoms[index]).leftNdx = index + 1;
        }

        GAPS_ASSERT_MSG(mAtoms[index].rightNdx != index + 1,
            index << " , " << size());
        GAPS_ASSERT_MSG(mAtoms[index].leftNdx != index + 1,
            index << " , " << size());
    }

    // delete atom from vector in O(1)
    mAtomPositions.erase(pos);
    mAtoms.pop_back();
    mUsedPositions.erase(pos);
    GAPS_ASSERT(mAtoms.size() == mAtomPositions.size());
    GAPS_ASSERT(mAtoms.size() == mUsedPositions.size());
}

// O(logN)
void AtomicDomain::updateMass(uint64_t pos, float newMass)
{
    GAPS_ASSERT_MSG(mAtomPositions.find(pos) != mAtomPositions.end(),
        pos);
    mAtoms[mAtomPositions.at(pos)].mass = newMass;
}

bool AtomicDomain::test(uint64_t pos) const
{
    return mAtomPositions.find(pos) != mAtomPositions.end();
}