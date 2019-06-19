#include "Atom.h"

Atom::Atom(uint64_t p, float m)
    : mLeft(NULL), mRight(NULL), mPos(p), mIterator(), mIndex(0), mMass(m)
{}

uint64_t Atom::pos() const
{
    return mPos;
}

float Atom::mass() const
{
    return mMass;
}

bool Atom::hasLeft() const
{
    return mLeft != NULL;
}

bool Atom::hasRight() const
{
    return mRight != NULL;
}

Atom* Atom::left() const
{
    return mLeft;
}

Atom* Atom::right() const
{
    return mRight;
}

void Atom::updateMass(float newMass)
{
    mMass = newMass;
}

void Atom::updatePos(uint64_t newPos)
{
    mPos = newPos;
}

void Atom::setLeft(Atom* atom)
{
    mLeft = atom;
}

void Atom::setRight(Atom* atom)
{
    mRight = atom;
}

void Atom::setIndex(unsigned index)
{
    mIndex = index;
}

void Atom::setIterator(AtomMapType::iterator it)
{
    mIterator = it;
}

unsigned Atom::index() const
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
