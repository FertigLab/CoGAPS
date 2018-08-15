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

void ProposalQueue::setDimensionSize(uint64_t binSize, uint64_t dimLength)
{
    mDimensionSize = binSize * dimLength;
    mUsedIndices.setDimensionSize(mNumBins / dimLength);
}

void ProposalQueue::setAlpha(float alpha)
{
    mAlpha = alpha;
}

void ProposalQueue::populate(AtomicDomain &domain, unsigned limit)
{
    unsigned nIter = 0;
    bool success = true;
    while (nIter++ < limit && success)
    {
        success = makeProposal(domain);
        if (!success)
        {
            mUseCachedRng = true;
        }
    }
}

void ProposalQueue::clear()
{
    mQueue.clear();
    mUsedPositions.clear();
    mUsedIndices.clear();
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);
}

unsigned ProposalQueue::size() const
{
    return mQueue.size();
}

const AtomicProposal& ProposalQueue::operator[](int n) const
{
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

void ProposalQueue::acceptBirth()
{
    #pragma omp atomic
    mMinAtoms++;
}

void ProposalQueue::rejectBirth()
{
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

    mU1 = mUseCachedRng ? mU1 : mRng.uniform();
    mU2 = mUseCachedRng ? mU2: mRng.uniform();
    mUseCachedRng = false;

    float lowerBound = deathProb(mMinAtoms);
    float upperBound = deathProb(mMaxAtoms);

    if (mU1 <= bdProb)
    {
        if (mU2 >= upperBound)
        {
            return birth(domain);
        }
        if (mU2 < lowerBound)
        {
            return death(domain);
        }
        return false;
    }
    return (mU1 < 0.75f || mMaxAtoms < 2) ? move(domain) : exchange(domain);
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

    if (!mUsedPositions.isEmptyInterval(lbound, rbound))
    {
        return false;
    }

    uint64_t newLocation = mRng.uniform64(lbound, rbound - 1);
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

    if (domain.hasRight(a1)) // has neighbor
    {
        if (!mUsedPositions.isEmptyInterval(a1.pos, a2.pos))
        {
            return false;
        }
    }
    else // exchange with first atom
    {
        if (!mUsedPositions.isEmptyInterval(a1.pos, mDomainSize))
        {
            return false;
        }
        
        if (!mUsedPositions.isEmptyInterval(0, domain.front().pos))
        {
            return false;
        }
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

Archive& operator<<(Archive &ar, ProposalQueue &queue)
{
    ar << queue.mMinAtoms << queue.mMaxAtoms << queue.mNumBins
        << queue.mDimensionSize << queue.mDomainSize << queue.mAlpha
        << queue.mRng;
    return ar;
}

Archive& operator>>(Archive &ar, ProposalQueue &queue)
{
    ar >> queue.mMinAtoms >> queue.mMaxAtoms >> queue.mNumBins
        >> queue.mDimensionSize >> queue.mDomainSize >> queue.mAlpha
        >> queue.mRng;
    return ar;
}
