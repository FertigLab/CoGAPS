#include "GapsAssert.h"
#include "ProposalQueue.h"
#include "Random.h"

void ProposalQueue::populate(const AtomicDomain &domain, unsigned limit)
{
    unsigned nIter = 0;
    while (nIter++ < limit && makeProposal(domain));
}

void ProposalQueue::setNumBins(unsigned nBins)
{
    mNumBins = nBins;
}

void ProposalQueue::setDomainSize(uint64_t size)
{
    mDomainSize = size;
}

void ProposalQueue::setDimensionSize(unsigned size)
{
    mDimensionSize = size;
}

void ProposalQueue::setAlpha(float alpha)
{
    mAlpha = alpha;
}

void ProposalQueue::clear(unsigned n)
{
    mQueue.erase(mQueue.end() - n, mQueue.end());
    mUsedIndices.erase(mUsedIndices.end() - n, mUsedIndices.end());
    mUsedPositions.erase(mUsedPositions.end() - n, mUsedPositions.end());
    GAPS_ASSERT(mMaxAtoms - mMinAtoms <= mQueue.size());
}

unsigned ProposalQueue::size() const
{
    return mQueue.size();
}

const AtomicProposal& ProposalQueue::operator[](int n) const
{
    return mQueue[mQueue.size() - 1 - n];
}

void ProposalQueue::acceptDeath()
{
    mMaxAtoms--;
}

void ProposalQueue::rejectDeath()
{
    mMinAtoms++;
}

void ProposalQueue::rejectBirth()
{
    mMinAtoms--;
    mMaxAtoms--;
}

float ProposalQueue::deleteProb(unsigned nAtoms) const
{
    double size = static_cast<double>(mDomainSize);
    double term1 = (size - static_cast<double>(nAtoms)) / size;
    double term2 = mAlpha * static_cast<double>(mNumBins) * term1;
    return static_cast<double>(nAtoms) / (static_cast<double>(nAtoms) + term2);
}

bool ProposalQueue::makeProposal(const AtomicDomain &domain)
{
    // special indeterminate cases
    if ((mMinAtoms == 0 && mMaxAtoms > 0) || (mMinAtoms < 2 && mMaxAtoms >= 2))
    {
        return false;
    }

    // always birth when no atoms exist
    if (mMaxAtoms == 0)
    {
        return birth(domain);
    }

    float bdProb = mMaxAtoms < 2 ? 0.6667f : 0.5f;
    float u1 = gaps::random::uniform(); // cache these values if a failure
    float u2 = gaps::random::uniform(); // happens, needed to prevent bias
    float lowerBound = deleteProb(mMinAtoms);
    float upperBound = deleteProb(mMaxAtoms);

    if (u1 < bdProb && u2 >= upperBound)
    {
        return birth(domain);
    }
    else if (u1 < bdProb && u2 < lowerBound)
    {
        return death(domain);
    }
    else if (u1 >= bdProb)
    {
        return (u1 < 0.75f || mMaxAtoms < 2) ? move(domain) : exchange(domain);
    }
    return false;
}

static bool isInVector(const std::vector<unsigned> &vec, unsigned n)
{
    return std::find(vec.begin(), vec.end(), n) != vec.end();
}

bool ProposalQueue::birth(const AtomicDomain &domain)
{
    uint64_t pos = domain.randomFreePosition();
    //if (isInVector(mUsedIndices, pos / mDimensionSize))
    //{
    //    return false; // matrix conflict - can't compute gibbs mass
    //}
    mQueue.push_back(AtomicProposal('B', pos));
    //mUsedIndices.push_back(pos / mDimensionSize);
    //mUsedPositions.push_back(pos);
    mMinAtoms++;
    mMaxAtoms++;
    return true;
}

bool ProposalQueue::death(const AtomicDomain &domain)
{
    Atom a = domain.randomAtom();
    //if (isInVector(mUsedIndices, a.pos / mDimensionSize))
    //{
    //    return false; // matrix conflict - can't compute gibbs mass or deltaLL
    //}
    mQueue.push_back(AtomicProposal('D', a.pos, a.mass));
    //mUsedIndices.push_back(a.pos / mDimensionSize);
    //mUsedPositions.push_back(a.pos);
    mMinAtoms--;
    return true;
}

bool ProposalQueue::move(const AtomicDomain &domain)
{
    Atom a = domain.randomAtom();
    uint64_t lbound = domain.hasLeft(a) ? domain.left(a).pos : 0;
    uint64_t rbound = domain.hasRight(a) ? domain.right(a).pos : mDomainSize;

    //for (unsigned i = 0; i < mUsedPositions.size(); ++i)
    //{
    //    if (mUsedPositions[i] >= lbound && mUsedPositions[i] <= rbound)
    //    {
    //        return false; // atomic conflict - don't know neighbors
    //    }
    //}

    uint64_t newLocation = gaps::random::uniform64(lbound + 1, rbound - 1);
    //if (isInVector(mUsedIndices, a.pos / mDimensionSize) || isInVector(mUsedIndices, newLocation / mDimensionSize))
    //{
    //    return false; // matrix conflict - can't compute deltaLL
    //}

    mQueue.push_back(AtomicProposal('M', a.pos, a.mass, newLocation));
    //mUsedIndices.push_back(a.pos / mDimensionSize);
    //mUsedIndices.push_back(newLocation / mDimensionSize);
    //mUsedPositions.push_back(a.pos);
    //mUsedPositions.push_back(newLocation);
    return true;
}

bool ProposalQueue::exchange(const AtomicDomain &domain)
{
    Atom a1 = domain.randomAtom();
    domain.test(a1.pos);
    GAPS_ASSERT(a1.rightNdx <= domain.size());
    GAPS_ASSERT(a1.leftNdx <= domain.size());
    //if (a1.right) // has neighbor
    //{
    //    for (unsigned i = 0; i < mUsedPositions.size(); ++i)
    //    {
    //        if (mUsedPositions[i] >= a1.pos && mUsedPositions[i] <= a1.right->pos)
    //        {
    //            return false; // atomic conflict - don't know right neighbor
    //        }
    //    }
    //}
    //else // exchange with first atom
    //{
    //    for (unsigned i = 0; i < mUsedPositions.size(); ++i)
    //    {
    //        if (mUsedPositions[i] >= a1.pos || mUsedPositions[i] <= domain.front().pos)
    //        {
    //            return false; // atomic conflict - don't know right neighbor
    //        }
    //    }
    //}

    Atom a2 = domain.hasRight(a1) ? domain.right(a1) : domain.front();
    domain.test(a1.pos);
    GAPS_ASSERT(a1.rightNdx <= domain.size());
    GAPS_ASSERT(a1.leftNdx <= domain.size());
    GAPS_ASSERT(a2.rightNdx <= domain.size());
    GAPS_ASSERT(a2.leftNdx <= domain.size());
    domain.test(a2.pos);
    //if (isInVector(mUsedIndices, a1.pos / mDimensionSize) || isInVector(mUsedIndices, a2.pos / mDimensionSize))
    //{
    //    return false; // matrix conflict - can't compute gibbs mass or deltaLL
    //}

    mQueue.push_back(AtomicProposal('E', a1.pos, a1.mass, a2.pos, a2.mass));
    //mUsedIndices.push_back(a1.pos / mDimensionSize);
    //mUsedIndices.push_back(a2.pos / mDimensionSize);
    //mUsedPositions.push_back(a1.pos);
    //mUsedPositions.push_back(a2.pos);
    mMinAtoms--;
    return true;
}