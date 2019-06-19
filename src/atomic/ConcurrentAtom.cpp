#include "ConcurrentAtom.h"

ConcurrentAtom::ConcurrentAtom(uint64_t p, float m)
    : mPos(p), mLeft(NULL), mRight(NULL), mIterator(), mIndex(0), mMass(m)
{}

uint64_t ConcurrentAtom::pos() const
{
    return mPos;
}

float ConcurrentAtom::mass() const
{
    return mMass;
}

bool ConcurrentAtom::hasLeft() const
{
    return mLeft != NULL;
}

bool ConcurrentAtom::hasRight() const
{
    return mRight != NULL;
}

ConcurrentAtom* ConcurrentAtom::left() const
{
    return mLeft;
}

ConcurrentAtom* ConcurrentAtom::right() const
{
    return mRight;
}

void ConcurrentAtom::updateMass(float newMass)
{
    mMass = newMass;
}

void ConcurrentAtom::updatePos(uint64_t newPos)
{
    mPos = newPos;
}

void ConcurrentAtom::setLeft(ConcurrentAtom* ConcurrentAtom)
{
    mLeft = ConcurrentAtom;
}

void ConcurrentAtom::setRight(ConcurrentAtom* ConcurrentAtom)
{
    mRight = ConcurrentAtom;
}

void ConcurrentAtom::setIndex(unsigned index)
{
    mIndex = index;
}

void ConcurrentAtom::setIterator(ConcurrentAtomMapType::iterator it)
{
    mIterator = it;
}

unsigned ConcurrentAtom::index() const
{
    return mIndex;
}

ConcurrentAtomMapType::iterator ConcurrentAtom::iterator() const
{
    return mIterator;
}

Archive& operator<<(Archive &ar, const ConcurrentAtom &a)
{
    ar << a.mPos << a.mMass;
    return ar;
}

Archive& operator>>(Archive &ar, ConcurrentAtom &a)
{
    ar >> a.mPos >> a.mMass;
    return ar;
}
