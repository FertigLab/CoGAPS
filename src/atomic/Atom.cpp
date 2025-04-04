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

Atom::Atom(uint64_t p, float m): mIterator(), mPos(p), mHasRight(false), mHasLeft(false), mMass(m)
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

void Atom::setLeftIndex(size_t index)
{
    mLeftIndex = index;
    mHasLeft = true;
}

void Atom::setRightIndex(size_t index)
{
    mRightIndex = index;
    mHasRight = true;
}

void Atom::setIndex(size_t index)
{
    mIndex = index;
}

void Atom::setIterator(AtomMapType::iterator it)
{
    mIterator = it;
}

bool Atom::hasLeft() const
{
    return mHasLeft;
}

bool Atom::hasRight() const
{
    return mHasRight;
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
