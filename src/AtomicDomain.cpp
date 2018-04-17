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
    iter = mAtomPositions.insert(std::pair<uint64_t, uint64_t>(pos, mAtoms.size())).first;
    iterLeft = iter;
    iterRight = iter;

    // find neighbors
    Atom atom(pos, mass);
    if (iter != mAtomPositions.begin())
    {
        --iterLeft;
        atom.leftNdx = iterLeft->second + 1;
    }
    if (++iter != mAtomPositions.end())
    {
        ++iterRight;
        atom.rightNdx = iterRight->second + 1;
    } 

    // add atom to vector
    GAPS_ASSERT(atom.rightNdx <= size());
    GAPS_ASSERT(atom.leftNdx <= size());
    mAtoms.push_back(atom);
    mUsedPositions.insert(pos);
    GAPS_ASSERT(mAtoms.size() == mAtomPositions.size());
    GAPS_ASSERT(mAtoms.size() == mUsedPositions.size());
}

void AtomicDomain::swap(uint64_t i1, uint64_t i2)
{
    if (i1 != i2)
    {
        // store the atom in position 1
        Atom temp1 = mAtoms[i1];
        Atom temp2 = mAtoms[i2];

        // switch atoms
        mAtoms[i1] = temp2;
        mAtoms[i2] = temp1;

        // update neighbors
        if (hasLeft(temp1))
        {
            GAPS_ASSERT(i2 + 1 <= size());
            left(temp1).rightNdx = i2 + 1;
        }
        if (hasRight(temp1))
        {
            GAPS_ASSERT(i2 + 1 <= size());
            right(temp1).leftNdx = i2 + 1;
        }
        if (hasLeft(temp2))
        {
            GAPS_ASSERT(i1 + 1 <= size());
            left(temp2).rightNdx = i1 + 1;
        }
        if (hasRight(temp2))
        {
            GAPS_ASSERT(i1 + 1 <= size());
            right(temp2).leftNdx = i1 + 1;
        }

        // update atom positions
        mAtomPositions.erase(mAtoms[i1].pos);
        mAtomPositions.erase(mAtoms[i2].pos);
        mAtomPositions.insert(std::pair<uint64_t, uint64_t>(mAtoms[i1].pos,
            i1));
        mAtomPositions.insert(std::pair<uint64_t, uint64_t>(mAtoms[i2].pos,
            i2));
    }
}

// O(logN)
// erase directly from position map
// swap with last atom in vector, pop off the back
void AtomicDomain::erase(uint64_t pos, bool display)
{
    // get vector index of this atom and erase it
    GAPS_ASSERT(mAtomPositions.find(pos) != mAtomPositions.end());
    uint64_t index = mAtomPositions.at(pos);

    // move atom to back
    swap(index, mAtoms.size() - 1);

    // connect neighbors of atom to be deleted
    if (hasLeft(mAtoms.back()))
    {
        left(mAtoms.back()).rightNdx = mAtoms.back().rightNdx;
    }
    if (hasRight(mAtoms.back()))
    {
        right(mAtoms.back()).leftNdx = mAtoms.back().leftNdx;
    }

    // delete atom from vector in O(1)
    GAPS_ASSERT(mAtomPositions.at(pos) == mAtoms.size() - 1);
    mAtomPositions.erase(pos);
    mAtoms.pop_back();
    mUsedPositions.erase(pos);
    GAPS_ASSERT(mAtoms.size() == mAtomPositions.size());
    GAPS_ASSERT(mAtoms.size() == mUsedPositions.size());
}

// O(logN)
void AtomicDomain::updateMass(uint64_t pos, float newMass)
{
    GAPS_ASSERT(mAtomPositions.find(pos) != mAtomPositions.end());
    mAtoms[mAtomPositions.at(pos)].mass = newMass;
}

void AtomicDomain::test(uint64_t pos) const
{
    GAPS_ASSERT(mAtomPositions.find(pos) != mAtomPositions.end());
}