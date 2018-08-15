#include "GapsAssert.h"
#include "AtomicDomain.h"
#include "math/Random.h"

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
    uint64_t ndx = mRng.uniform64(0, mAtoms.size() - 1);
    return mAtoms[ndx];
}

// Average Case O(1)
uint64_t AtomicDomain::randomFreePosition() const
{
    uint64_t pos = 0;
    do
    {
        pos = mRng.uniform64(0, mDomainSize);
    } while (mPositionLookup.count(pos) > 0);
    return pos;
}

// O(1)
uint64_t AtomicDomain::size() const
{
    return mAtoms.size();
}

// O(1)
Atom& AtomicDomain::_left(const Atom &atom)
{
    return mAtoms[atom.leftNdx - 1];
}

// O(1)
Atom& AtomicDomain::_right(const Atom &atom)
{
    return mAtoms[atom.rightNdx - 1];
}

// O(1)
const Atom& AtomicDomain::left(const Atom &atom) const
{
    return mAtoms[atom.leftNdx - 1];
}

// O(1)
const Atom& AtomicDomain::right(const Atom &atom) const
{
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
Atom AtomicDomain::insert(uint64_t pos, float mass)
{
    // insert position into map
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
        _left(atom).rightNdx = size() + 1;
    }
    if (++iter != mAtomPositions.end())
    {
        ++iterRight;
        atom.rightNdx = iterRight->second + 1;
        _right(atom).leftNdx = size() + 1;
    } 

    // add atom to vector and lookup table
    mPositionLookup.insert(std::pair<uint64_t, uint64_t>(pos, size()));
    mAtoms.push_back(atom);

    return atom;
}

// O(logN)
// erase directly from position map and used positions hash set
// swap with last atom in vector, pop off the back
void AtomicDomain::erase(uint64_t pos)
{
    // get vector index of this atom
    uint64_t index = mPositionLookup.at(pos);

    // connect neighbors of atom to be deleted
    if (hasLeft(mAtoms[index]))
    {
        _left(mAtoms[index]).rightNdx = mAtoms[index].rightNdx;
    }
    if (hasRight(mAtoms[index]))
    {
        _right(mAtoms[index]).leftNdx = mAtoms[index].leftNdx;
    }

    // replace with atom from back
    if (index < size() - 1)
    {
        mAtoms[index] = mAtoms.back();

        // update position in map
        mAtomPositions.erase(mAtoms[index].pos);
        mAtomPositions.insert(std::pair<uint64_t, uint64_t>(mAtoms[index].pos,
            index));

        mPositionLookup.erase(mAtoms[index].pos);
        mPositionLookup.insert(std::pair<uint64_t, uint64_t>(mAtoms[index].pos,
            index));
    
        // update moved atom's neighbors
        if (hasLeft(mAtoms[index]))
        {
            _left(mAtoms[index]).rightNdx = index + 1;
        }
        if (hasRight(mAtoms[index]))
        {
            _right(mAtoms[index]).leftNdx = index + 1;
        }
    }

    // delete atom from vector in O(1), map in O(logN)
    mAtomPositions.erase(pos);
    mAtoms.pop_back();
    mPositionLookup.erase(pos);
}

void AtomicDomain::cacheInsert(uint64_t pos, float mass) const
{
    unsigned ndx = 0;
    #pragma omp critical(atomicInsert)
    {
        ndx = mInsertCacheIndex++;
    }
    mInsertCache[ndx] = RawAtom(pos, mass);
}

void AtomicDomain::cacheErase(uint64_t pos) const
{
    unsigned ndx = 0;
    #pragma omp critical(atomicErase)
    {
        ndx = mEraseCacheIndex++;
    }
    mEraseCache[ndx] = pos;
}

void AtomicDomain::resetCache(unsigned n)
{
    mInsertCacheIndex = 0;
    mEraseCacheIndex = 0;
    mInsertCache.resize(n);
    mEraseCache.resize(n);
}

void AtomicDomain::flushCache()
{
    for (unsigned i = 0; i < mEraseCacheIndex; ++i)
    {
        erase(mEraseCache[i]);
    }

    for (unsigned i = 0; i < mInsertCacheIndex; ++i)
    {
        insert(mInsertCache[i].pos, mInsertCache[i].mass);
    }

    mInsertCache.clear();
    mEraseCache.clear();
}

// O(1)
void AtomicDomain::updateMass(uint64_t pos, float newMass)
{
    mAtoms[mPositionLookup.at(pos)].mass = newMass;
}

Archive& operator<<(Archive &ar, Atom &a)
{
    ar << a.leftNdx << a.rightNdx << a.pos << a.mass;
    return ar;
}

Archive& operator>>(Archive &ar, Atom &a)
{
    ar >> a.leftNdx >> a.rightNdx >> a.pos >> a.mass;
    return ar;
}

Archive& operator<<(Archive &ar, AtomicDomain &domain)
{
    unsigned nAtoms = domain.mAtoms.size();
    ar << nAtoms << domain.mDomainSize;

    for (unsigned i = 0; i < nAtoms; ++i)
    {
        ar << domain.mAtoms[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, AtomicDomain &domain)
{
    unsigned nAtoms = 0;
    ar >> nAtoms >> domain.mDomainSize;

    Atom a(0, 0.f);
    for (unsigned i = 0; i < nAtoms; ++i)
    {
        ar >> a;
        domain.insert(a.pos, a.mass);
    }
    return ar;
}
