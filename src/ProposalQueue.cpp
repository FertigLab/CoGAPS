#include "GapsAssert.h"
#include "ProposalQueue.h"
#include "math/Random.h"

//////////////////////////////// AtomicProposal ////////////////////////////////

// birth
AtomicProposal::AtomicProposal(char t, uint64_t pos)
    : type(t), birthPos(pos), moveDest(0), atom1(NULL), atom2(NULL)
{}
    
// death
AtomicProposal::AtomicProposal(char t, Atom *atom)
    : type(t), birthPos(0), moveDest(0), atom1(atom), atom2(NULL)
{}

// move
AtomicProposal::AtomicProposal(char t, Atom *atom, uint64_t dest)
    : type(t), birthPos(0), moveDest(dest), atom1(atom), atom2(NULL)
{}

// exchange
AtomicProposal::AtomicProposal(char t, Atom *a1, Atom *a2)
    : type(t), birthPos(0), moveDest(0), atom1(a1), atom2(a2)
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

float ProposalQueue::deathProb(uint64_t nAtoms) const
{
    double size = static_cast<double>(mDomainLength);
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
        return false; // can't determine B/D since range is too wide
    }
    return (mU1 < 0.75f || mMaxAtoms < 2) ? move(domain) : exchange(domain);
}
    
unsigned ProposalQueue::primaryIndex(uint64_t pos) const
{
    return pos / mSecondaryDimLength;
}

unsigned ProposalQueue::secondaryIndex(uint64_t pos) const
{
    return (pos / mBinLength) % mSecondaryDimSize;
}

// TODO add atoms with empty mass? fill in mass in gibbssampler?
// inserting invalidates previous pointers, but not inserting
// prevents them from being selected for death
bool ProposalQueue::birth(AtomicDomain &domain)
{
    uint64_t pos = domain.randomFreePosition();
    if (mUsedIndices.contains(primaryIndex(pos)))
    {
        return false; // matrix conflict - can't compute gibbs mass
    }

    mQueue.push_back(AtomicProposal('B', pos));
    mUsedIndices.insert(primaryIndex(pos));
    mUsedPositions.insert(pos);
    ++mMaxAtoms;
    return true;
}

bool ProposalQueue::death(AtomicDomain &domain)
{
    Atom* a = domain.randomAtom();
    if (mUsedIndices.contains(primaryIndex(a->pos)))
    {
        return false; // matrix conflict - can't compute gibbs mass or deltaLL
    }

    mQueue.push_back(AtomicProposal('D', a));
    mUsedIndices.insert(primaryIndex(a->pos));
    mUsedPositions.insert(a->pos);
    --mMinAtoms;
    return true;
}

bool ProposalQueue::move(AtomicDomain &domain)
{
    AtomNeighborhood hood = domain.randomAtomWithNeighbors();
    uint64_t lbound = hood.hasLeft() ? hood.left->pos : 0;
    uint64_t rbound = hood.hasRight() ? hood.right->pos : mDomainLength;

    if (!mUsedPositions.isEmptyInterval(lbound, rbound))
    {
        return false; // atomic conflict - don't know neighbors
    }

    uint64_t newLocation = mRng.uniform64(lbound + 1, rbound - 1);

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
        return true; // TODO automatically accept exchanges in same bin
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
