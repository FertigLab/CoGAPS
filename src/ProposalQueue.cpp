#include "GapsAssert.h"
#include "ProposalQueue.h"
#include "math/Random.h"

void ProposalQueue::setNumBins(unsigned nBins)
{
    mNumBins = nBins;
}

void ProposalQueue::setDomainSize(uint64_t size)
{
    mDomainSize = size;
}

void ProposalQueue::setDimensionSize(uint64_t size)
{
    mDimensionSize = size;
}

void ProposalQueue::setAlpha(float alpha)
{
    mAlpha = alpha;
}

void ProposalQueue::populate(AtomicDomain &domain, unsigned limit)
{
    unsigned nIter = 0;
    while (nIter++ < limit && makeProposal(domain));
}

void ProposalQueue::clear(unsigned n)
{
    //mQueue.erase(mQueue.end() - n, mQueue.end());
    //mUsedIndices.erase(mUsedIndices.end() - n, mUsedIndices.end());
    //mUsedPositions.erase(mUsedPositions.end() - n, mUsedPositions.end());
    mQueue.clear();
    mUsedPositions.clear();
    mUsedIndices.clear();
    //GAPS_ASSERT(mMaxAtoms - mMinAtoms <= mQueue.size());
}

unsigned ProposalQueue::size() const
{
    return mQueue.size();
}

const AtomicProposal& ProposalQueue::operator[](int n) const
{
    //return mQueue[mQueue.size() - 1 - n];
    return mQueue[n];
}

void ProposalQueue::acceptDeath()
{
    #pragma omp atomic
    mMaxAtoms--;
}

void ProposalQueue::rejectDeath()
{
    #pragma omp atomic
    mMinAtoms++;
}

void ProposalQueue::rejectBirth()
{
    #pragma omp atomic
    mMinAtoms--;

    #pragma omp atomic
    mMaxAtoms--;
}

float ProposalQueue::deathProb(uint64_t nAtoms) const
{
    double size = static_cast<double>(mDomainSize);
    double term1 = (size - static_cast<double>(nAtoms)) / size;
    double term2 = mAlpha * static_cast<double>(mNumBins) * term1;
    return static_cast<double>(nAtoms) / (static_cast<double>(nAtoms) + term2);
}

bool ProposalQueue::makeProposal(AtomicDomain &domain)
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
    float lowerBound = deathProb(mMinAtoms);
    float upperBound = deathProb(mMaxAtoms);

    if (u1 <= bdProb)
    {
        float u2 = gaps::random::uniform(); // happens, needed to prevent bias
        if (u2 >= upperBound)
        {
            return birth(domain);
        }
        if (u2 < lowerBound)
        {
            return death(domain);
        }
        return false;
    }
    else if (u1 >= bdProb)
    {
        return (u1 < 0.75f || mMaxAtoms < 2) ? move(domain) : exchange(domain);
    }
    return false;
}

static bool isInVector(const std::vector<uint64_t> &vec, uint64_t n)
{
    return std::find(vec.begin(), vec.end(), n) != vec.end();
}

bool ProposalQueue::birth(AtomicDomain &domain)
{
    uint64_t pos = domain.randomFreePosition();
    if (mUsedIndices.count(pos / mDimensionSize))
    {
        return false; // matrix conflict - can't compute gibbs mass
    }

    mQueue.push_back(AtomicProposal('B', pos));
    mUsedIndices.insert(pos / mDimensionSize);
    mUsedPositions.insert(pos);
    domain.insert(pos, 0.f);
    ++mMinAtoms;
    ++mMaxAtoms;
    return true;
}

bool ProposalQueue::death(AtomicDomain &domain)
{
    Atom a = domain.randomAtom();
    if (mUsedIndices.count(a.pos / mDimensionSize))
    {
        return false; // matrix conflict - can't compute gibbs mass or deltaLL
    }

    mQueue.push_back(AtomicProposal('D', a.pos, a.mass));
    mUsedIndices.insert(a.pos / mDimensionSize);
    mUsedPositions.insert(a.pos);
    --mMinAtoms;
    return true;
}

bool ProposalQueue::move(AtomicDomain &domain)
{
    Atom a = domain.randomAtom();
    uint64_t lbound = domain.hasLeft(a) ? domain.left(a).pos : 0;
    uint64_t rbound = domain.hasRight(a) ? domain.right(a).pos : mDomainSize;

    uint64_t low = *(mUsedPositions.lower_bound(lbound));
    if (low != *(mUsedPositions.end()) && low <= rbound)
    {
        return false; //atomic conflict - don't know neighbors
    }

    uint64_t newLocation = gaps::random::uniform64(lbound, rbound - 1);
    if (mUsedIndices.count(a.pos / mDimensionSize) || mUsedIndices.count(newLocation / mDimensionSize))
    {
        return false; // matrix conflict - can't compute deltaLL
    }

    mQueue.push_back(AtomicProposal('M', a.pos, a.mass, newLocation));
    mUsedIndices.insert(a.pos / mDimensionSize);
    mUsedIndices.insert(newLocation / mDimensionSize);
    mUsedPositions.insert(a.pos);
    mUsedPositions.insert(newLocation);
    return true;
}

bool ProposalQueue::exchange(AtomicDomain &domain)
{
    Atom a1 = domain.randomAtom();
    Atom a2 = domain.hasRight(a1) ? domain.right(a1) : domain.front();

    //if (domain.hasRight(a1))
    //{
    //    if (*(mUsedPositions.lower_bound(a1.pos)) <= a2.pos)
    //    {
    //        return false;
    //    }
    //}
    //else
    //{
    //    if (*(mUsedPositions.end()) >= a1.pos || *(mUsedPositions.begin()) <= a2.pos)
    //    {
    //        return false;
    //    }
    //}

    if (domain.hasRight(a1)) // has neighbor
    {
        std::set<uint64_t>::iterator low = mUsedPositions.lower_bound(a1.pos);
        if (low != mUsedPositions.end() && *low <= a2.pos)
        {
            return false;
        }
    }
    else // exchange with first atom
    {
        for (std::set<uint64_t>::iterator it = mUsedPositions.begin(); it != mUsedPositions.end(); ++it)
        {
            if (*it >= a1.pos || *it <= domain.front().pos)
            {
                return false; // atomic conflict - don't know right neighbor
            }
        }
        //std::set<uint64_t>::iterator it = mUsedPositions.upper_bound(a1.pos);
        //if (it != mUsedPositions.end() || *(mUsedPositions.begin()) <= a2.pos)
        //{
        //    return false;
        //}
    }

    if (mUsedIndices.count(a1.pos / mDimensionSize) || mUsedIndices.count(a2.pos / mDimensionSize))
    {
        return false; // matrix conflict - can't compute gibbs mass or deltaLL
    }

    mQueue.push_back(AtomicProposal('E', a1.pos, a1.mass, a2.pos, a2.mass));
    mUsedIndices.insert(a1.pos / mDimensionSize);
    mUsedIndices.insert(a2.pos / mDimensionSize);
    mUsedPositions.insert(a1.pos);
    mUsedPositions.insert(a2.pos);
    --mMinAtoms;
    return true;
}
