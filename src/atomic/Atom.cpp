#include "Atom.h"
#include "../utils/Archive.h"

#include <cstddef>

AtomNeighborhood::AtomNeighborhood()
    : center(NULL), left(NULL), right(NULL)
{}

AtomNeighborhood::AtomNeighborhood(Atom *l, Atom *c, Atom *r)
    : center(c), left(l), right(r)
{}

bool AtomNeighborhood::hasLeft() const
{
    return left != NULL;
}

bool AtomNeighborhood::hasRight() const
{
    return right != NULL;
}

Atom::Atom(uint64_t p, float m)
    : mIterator(), mPos(p), mLeftIndex(-1), mRightIndex(-1), mIndex(-1), mMass(m)
{}

uint64_t Atom::pos() const
{
    return mPos;
}

float Atom::mass() const
{
    return mMass;
}

void Atom::updateMass(float newMass)
{
    mMass = newMass;
}

void Atom::updatePos(uint64_t newPos)
{
    mPos = newPos;
}

void Atom::setLeftIndex(int index)
{
    mLeftIndex = index;
}

void Atom::setRightIndex(int index)
{
    mRightIndex = index;
}

void Atom::setIndex(int index)
{
    mIndex = index;
}

void Atom::setIterator(AtomMapType::iterator it)
{
    mIterator = it;
}

bool Atom::hasLeft() const
{
    return mLeftIndex >= 0;
}

bool Atom::hasRight() const
{
    return mRightIndex >= 0;
}

int Atom::leftIndex() const
{
    return mLeftIndex;
}

int Atom::rightIndex() const
{
    return mRightIndex;
}

int Atom::index() const
{
    return mIndex;
}

AtomMapType::iterator Atom::iterator() const
{
    return mIterator;
}

Archive& operator<<(Archive &ar, const Atom &a)
{
    ar << a.mPos << a.mMass;
    return ar;
}

Archive& operator>>(Archive &ar, Atom &a)
{
    ar >> a.mPos >> a.mMass;
    return ar;
}
