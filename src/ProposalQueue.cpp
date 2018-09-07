#include "GapsAssert.h"
#include "ProposalQueue.h"
#include "math/Random.h"

//////////////////////////////// AtomicProposal ////////////////////////////////

// birth/death
AtomicProposal::AtomicProposal(char t)
    : pos(0), atom1(NULL), atom2(NULL), type(t)
{}
    
//////////////////////////////// ProposalQueue /////////////////////////////////

ProposalQueue::ProposalQueue(unsigned primaryDimSize, unsigned secondaryDimSize)
    :
mMinAtoms(0), mMaxAtoms(0), mNumBins(primaryDimSize * secondaryDimSize),
mBinLength(std::numeric_limits<uint64_t>::max() / mNumBins),
mSecondaryDimLength(mBinLength * secondaryDimSize),
mDomainLength(mBinLength * mNumBins), mSecondaryDimSize(secondaryDimSize),
mAlpha(0.f), mU1(0.f), mU2(0.f), mUseCachedRng(false)
{
    mUsedIndices.setDimensionSize(primaryDimSize);
}

void ProposalQueue::setAlpha(float alpha)
{
    mAlpha = alpha;
}

void ProposalQueue::populate(AtomicDomain &domain, unsigned limit)
{
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);
    GAPS_ASSERT(mMaxAtoms == domain.size());

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
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);

    mQueue.clear();
    mUsedPositions.clear();
    mUsedIndices.clear();
}

unsigned ProposalQueue::size() const
{
    return mQueue.size();
}

AtomicProposal& ProposalQueue::operator[](int n)
{
    GAPS_ASSERT(mQueue.size() > 0);
    GAPS_ASSERT(n < mQueue.size());

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

float ProposalQueue::deathProb(double nAtoms) const
{
    double numer = nAtoms * mDomainLength;
    return numer / (numer + mAlpha * mNumBins * (mDomainLength - nAtoms));
}

bool ProposalQueue::makeProposal(AtomicDomain &domain)
{
    // special indeterminate case
    if (mMinAtoms < 2 && mMaxAtoms >= 2)
    {
        return false;
    }

    // always birth when no atoms exist
    if (mMaxAtoms < 2)
    {
        return birth(domain);
    }

    mU1 = mUseCachedRng ? mU1 : mRng.uniform();
    mU2 = mUseCachedRng ? mU2: mRng.uniform();
    mUseCachedRng = false;

    float lowerBound = deathProb(static_cast<double>(mMinAtoms));
    float upperBound = deathProb(static_cast<double>(mMaxAtoms));

    if (mU1 < 0.5f)
    {
        if (mU2 < lowerBound)
        {
            return death(domain);
        }
        if (mU2 >= upperBound)
        {
            return birth(domain);
        }
        return false; // can't determine B/D since range is too wide
    }
    return (mU1 < 0.75f) ? move(domain) : exchange(domain);
}
    
bool ProposalQueue::birth(AtomicDomain &domain)
{
    AtomicProposal prop('B');
    prop.r1 = (prop.atom1->pos / mBinSize) / mNumCols;
    prop.c1 = (prop.atom1->pos / mBinSize) % mNumCols;

    prop.init(); // initialize rng

    if (mUsedMatrixIndices.contains(prop.r1))
    {
        return false; // matrix conflict - can't compute alpha parameters
    }
    prop.atom1 = mDomain.insert(domain.randomFreePosition(), 0.f);

    mUsedMatrixIndices.insert(prop.r1);
    mUsedAtoms.insert(prop.atom1->pos);

    mQueue.push_back(prop);
    ++mMaxAtoms;
    return true;
}

bool ProposalQueue::death(AtomicDomain &domain)
{
    AtomicProposal prop('D');
    prop->atom1 = domain.randomAtom();
    prop.r1 = (prop.atom1->pos / mBinSize) / mNumCols;
    prop.c1 = (prop.atom1->pos / mBinSize) % mNumCols;

    if (mUsedMatrixIndices.contains(prop.r1))
    {
        return false; // matrix conflict - can't compute alpha parameters
    }
    mUsedMatrixIndices.insert(prop.r1);
    mUsedAtoms.insert(prop.atom1->pos);

    prop.init(); // initialize rng
    mQueue.push_back(prop);
    --mMinAtoms;
    return true;
}

bool ProposalQueue::move(AtomicDomain &domain)
{
    AtomicProposal prop('M');

    AtomNeighborhood hood = domain.randomAtomWithNeighbors();
    
    uint64_t lbound = hood.hasLeft() ? hood.left->pos : 0;
    uint64_t rbound = hood.hasRight() ? hood.right->pos : mDomainLength;

    uint64_t newLocation = mRng.uniform64(lbound + 1, rbound - 1);

    if (mUsedAtoms.contains(lbound) || mUsedAtoms.contains(rbound))
    {
        return false; // atomic conflict - don't know neighbors
    }

    if (primaryIndex(hood.center->pos) == primaryIndex(newLocation)
    && secondaryIndex(hood.center->pos) == secondaryIndex(newLocation))
    {
        hood.center->pos = newLocation; // automatically accept moves in same bin
        return true;
    }

    if (mUsedIndices.contains(primaryIndex(hood.center->pos))
    || mUsedIndices.contains(primaryIndex(newLocation)))
    {
        return false; // matrix conflict - can't compute deltaLL
    }

    mQueue.push_back(AtomicProposal('M', hood.center, newLocation));
    mUsedIndices.insert(primaryIndex(hood.center->pos));
    mUsedIndices.insert(primaryIndex(newLocation));
    mUsedPositions.insert(hood.center->pos);
    mUsedPositions.insert(newLocation);
    return true;
}

bool ProposalQueue::exchange(AtomicDomain &domain)
{
    AtomNeighborhood hood = domain.randomAtomWithRightNeighbor();
    Atom* a1 = hood.center;
    Atom* a2 = hood.hasRight() ? hood.right : domain.front();

    if (hood.hasRight()) // has neighbor
    {
        if (!mUsedPositions.isEmptyInterval(a1->pos, a2->pos))
        {
            return false; // atomic conflict - don't know right neighbor
        }
    }
    else // exchange with first atom
    {
        if (!mUsedPositions.isEmptyInterval(a1->pos, mDomainLength))
        {
            return false; // atomic conflict - don't know right neighbor
        }
        
        if (!mUsedPositions.isEmptyInterval(0, domain.front()->pos))
        {
            return false; // atomic conflict - don't know right neighbor
        }
    }

    if (primaryIndex(a1->pos) == primaryIndex(a2->pos)
    && secondaryIndex(a1->pos) == secondaryIndex(a2->pos))
    {
        GapsRng rng;
        float newMass = rng.truncGammaUpper(a1->mass + a2->mass, 2.f, 1.f / mLambda);
        float delta = (a1->mass > a2->mass) ? newMass - a1->mass : a2->mass - newMass;
        if (a1->mass + delta > gaps::epsilon && a2->mass - delta > gaps::epsilon)
        {
            a1->mass += delta;
            a2->mass -= delta;
        }        
        return true;
    }

    if (mUsedIndices.contains(primaryIndex(a1->pos))
    || mUsedIndices.contains(primaryIndex(a2->pos)))
    {
        return false; // matrix conflict - can't compute gibbs mass or deltaLL
    }

    mQueue.push_back(AtomicProposal('E', a1, a2));
    mUsedIndices.insert(primaryIndex(a1->pos));
    mUsedIndices.insert(primaryIndex(a2->pos));
    mUsedPositions.insert(a1->pos);
    mUsedPositions.insert(a2->pos);
    --mMinAtoms;
    return true;
}

Archive& operator<<(Archive &ar, ProposalQueue &q)
{
    ar << q.mMinAtoms << q.mMaxAtoms << q.mNumBins << q.mBinLength
        << q.mSecondaryDimLength << q.mDomainLength << q.mSecondaryDimSize
        << q.mAlpha << q.mRng;
    return ar;
}

Archive& operator>>(Archive &ar, ProposalQueue &q)
{
    ar >> q.mMinAtoms >> q.mMaxAtoms >> q.mNumBins >> q.mBinLength
        >> q.mSecondaryDimLength >> q.mDomainLength >> q.mSecondaryDimSize
        >> q.mAlpha >> q.mRng;
    return ar;
}