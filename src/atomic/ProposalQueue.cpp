#include "ProposalQueue.h"
#include "../atomic/ConcurrentAtomicDomain.h"
#include "../math/Math.h"
#include "../math/Random.h"
#include "../utils/Archive.h"
#include "../utils/GapsAssert.h"

#include <limits>

//////////////////////////////// AtomicProposal ////////////////////////////////

AtomicProposal::AtomicProposal(char t, GapsRandomState *randState)
    : rng(randState), pos(0), atom1(NULL), atom2(NULL), r1(0), c1(0), r2(0),
    c2(0), type(t)
{}
    
//////////////////////////////// ProposalQueue /////////////////////////////////

ProposalQueue::ProposalQueue(uint64_t nElements, uint64_t nPatterns,
GapsRandomState *randState)
    :
mUsedMatrixIndices(nElements / nPatterns),
mRandState(randState),
mRng(randState),
mMinAtoms(0),
mMaxAtoms(0),
mBinLength(std::numeric_limits<uint64_t>::max() / nElements),
mNumCols(nPatterns),
mAlpha(0.0),
mDomainLength(static_cast<double>(mBinLength * nElements)),
mNumBins(static_cast<double>(nElements)),
mU1(0.f),
mU2(0.f),
mNumProcessed(0),
mUseCachedRng(false)
{}

void ProposalQueue::setAlpha(float alpha)
{
    mAlpha = static_cast<double>(alpha);
}

void ProposalQueue::setLambda(float lambda)
{
    mLambda = lambda;
}

unsigned ProposalQueue::nProcessed() const
{
    return mNumProcessed;
}

void ProposalQueue::populate(ConcurrentAtomicDomain &domain, unsigned limit)
{
    GAPS_ASSERT(mQueue.empty());
    GAPS_ASSERT(mUsedAtoms.isEmpty());
    GAPS_ASSERT(mUsedMatrixIndices.isEmpty());
    GAPS_ASSERT(mProposedMoves.isEmpty());
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);
    GAPS_ASSERT_MSG(mMaxAtoms == domain.size(), mMaxAtoms << " != " << domain.size());

    bool success = true;
    mNumProcessed = 0;
    while (mNumProcessed < limit && success)
    {
        if (!makeProposal(domain))
        {
            success = false;
            mUseCachedRng = true;
        }
        else
        {
            ++mNumProcessed;
        }
    }
}

void ProposalQueue::clear()
{
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);
    mQueue.clear();
    mUsedMatrixIndices.clear();
    mUsedAtoms.clear();
    mProposedMoves.clear();
}

unsigned ProposalQueue::size() const
{
    return mQueue.size();
}

AtomicProposal& ProposalQueue::operator[](int n)
{
    GAPS_ASSERT(mQueue.size() > 0);
    GAPS_ASSERT(static_cast<unsigned>(n) < mQueue.size());
    return mQueue[n];
}

void ProposalQueue::acceptDeath()
{
    #pragma omp atomic
    --mMaxAtoms;
}

void ProposalQueue::rejectDeath()
{
    #pragma omp atomic
    ++mMinAtoms;
}

void ProposalQueue::acceptBirth()
{
    #pragma omp atomic
    ++mMinAtoms;
}

void ProposalQueue::rejectBirth()
{
    #pragma omp atomic
    --mMaxAtoms;
}

float ProposalQueue::deathProb(double nAtoms) const
{
    double numer = nAtoms * mDomainLength;
    return numer / (numer + mAlpha * mNumBins * (mDomainLength - nAtoms));
}

bool ProposalQueue::makeProposal(ConcurrentAtomicDomain &domain)
{
    mU1 = mUseCachedRng ? mU1 : mRng.uniform();
    mU2 = mUseCachedRng ? mU2: mRng.uniform();
    mUseCachedRng = false;

    if (mMinAtoms < 2 && mMaxAtoms >= 2)
    {
        return false; // special indeterminate case
    }

    if (mMaxAtoms < 2)
    {
        return birth(domain); // always birth when 0 or 1 atoms exist
    }

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

bool ProposalQueue::birth(ConcurrentAtomicDomain &domain)
{
    AtomicProposal prop('B', mRandState);
    uint64_t pos = domain.randomFreePosition(&(prop.rng));

    if (mProposedMoves.overlap(pos))
    {
        mRandState->rollBackOnce(); // ensure same proposal next time
        return false; // this birth would break assumption moves doesn't re-order domain
    }

    prop.r1 = (pos / mBinLength) / mNumCols;
    prop.c1 = (pos / mBinLength) % mNumCols;
    if (mUsedMatrixIndices.contains(prop.r1))
    {
        mRandState->rollBackOnce(); // ensure same proposal next time
        return false; // matrix conflict - can't compute gibbs mass
    }
    prop.atom1 = domain.insert(pos, 0.f);

    mUsedMatrixIndices.insert(prop.r1);
    mUsedAtoms.insert(prop.atom1->pos());
    mQueue.push_back(prop);
    ++mMaxAtoms;
    return true;
}

bool ProposalQueue::death(ConcurrentAtomicDomain &domain)
{
    AtomicProposal prop('D', mRandState);
    prop.atom1 = domain.randomAtom(&(prop.rng));
    prop.r1 = (prop.atom1->pos() / mBinLength) / mNumCols;
    prop.c1 = (prop.atom1->pos() / mBinLength) % mNumCols;

    if (mUsedMatrixIndices.contains(prop.r1))
    {
        mRandState->rollBackOnce(); // ensure same proposal next time
        return false; // matrix conflict - can't compute gibbs mass or deltaLL
    }

    mUsedMatrixIndices.insert(prop.r1);
    mUsedAtoms.insert(prop.atom1->pos());
    mQueue.push_back(prop);
    --mMinAtoms;
    return true;
}

bool ProposalQueue::move(ConcurrentAtomicDomain &domain)
{
    AtomicProposal prop('M', mRandState);
    ConcurrentAtomNeighborhood hood = domain.randomAtomWithNeighbors(&(prop.rng));
    prop.atom1 = hood.center;

    uint64_t lbound = hood.hasLeft() ? hood.left->pos() : 0;
    uint64_t rbound = hood.hasRight() ? hood.right->pos() : static_cast<uint64_t>(mDomainLength);

    if (mUsedAtoms.contains(lbound) || mUsedAtoms.contains(rbound))
    {
        mRandState->rollBackOnce(); // ensure same proposal next time
        return false; // atomic conflict - don't know neighbors
    }

    prop.pos = prop.rng.uniform64(lbound + 1, rbound - 1);
    prop.r1 = (prop.atom1->pos() / mBinLength) / mNumCols;
    prop.c1 = (prop.atom1->pos() / mBinLength) % mNumCols;
    prop.r2 = (prop.pos / mBinLength) / mNumCols;
    prop.c2 = (prop.pos / mBinLength) % mNumCols;

    if (mUsedMatrixIndices.contains(prop.r1) || mUsedMatrixIndices.contains(prop.r2))
    {
        mRandState->rollBackOnce(); // ensure same proposal next time
        return false; // matrix conflict - can't compute deltaLL
    }

    if (prop.r1 == prop.r2 && prop.c1 == prop.c2)
    {
        domain.move(prop.atom1, prop.pos);
        return true; // automatically accept moves in same bin
    }

    mQueue.push_back(prop);
    mUsedMatrixIndices.insert(prop.r1);
    mUsedMatrixIndices.insert(prop.r2);
    mUsedAtoms.insert(prop.atom1->pos());
    mProposedMoves.insert(prop.atom1->pos(), prop.pos);
    return true;
}

bool ProposalQueue::exchange(ConcurrentAtomicDomain &domain)
{
    AtomicProposal prop('E', mRandState);
    ConcurrentAtomNeighborhood hood = domain.randomAtomWithNeighbors(&(prop.rng));
    prop.atom1 = hood.center;
    prop.atom2 = hood.hasRight() ? hood.right : domain.front();
    prop.r1 = (prop.atom1->pos() / mBinLength) / mNumCols;
    prop.c1 = (prop.atom1->pos() / mBinLength) % mNumCols;
    prop.r2 = (prop.atom2->pos() / mBinLength) / mNumCols;
    prop.c2 = (prop.atom2->pos() / mBinLength) % mNumCols;

    if (mUsedMatrixIndices.contains(prop.r1) || mUsedMatrixIndices.contains(prop.r2))
    {
        mRandState->rollBackOnce(); // ensure same proposal next time
        return false; // matrix conflict - can't compute deltaLL or gibbs mass
    }

    if (prop.r1 == prop.r2 && prop.c1 == prop.c2)
    {
        float newMass = prop.rng.truncGammaUpper(prop.atom1->mass() + prop.atom2->mass(), 1.f / mLambda);
        float delta = (prop.atom1->mass() > prop.atom2->mass()) ? newMass - prop.atom1->mass() : prop.atom2->mass() - newMass;
        if (prop.atom1->mass() + delta > gaps::epsilon && prop.atom2->mass() - delta > gaps::epsilon)
        {
            prop.atom1->updateMass(prop.atom1->mass() + delta);
            prop.atom2->updateMass(prop.atom2->mass() - delta);
        }        
        return true; // automatically accept exchanges in same bin
    }

    mQueue.push_back(prop);
    mUsedMatrixIndices.insert(prop.r1);
    mUsedMatrixIndices.insert(prop.r2);
    return true;
}

Archive& operator<<(Archive &ar, const ProposalQueue &q)
{
    ar << q.mRng << q.mMinAtoms << q.mMaxAtoms << q.mBinLength << q.mNumCols
        << q.mAlpha << q.mDomainLength << q.mNumBins << q.mLambda
        << q.mUseCachedRng << q.mU1 << q.mU2;
    return ar;
}

Archive& operator>>(Archive &ar, ProposalQueue &q)
{
    ar >> q.mRng >> q.mMinAtoms >> q.mMaxAtoms >> q.mBinLength >> q.mNumCols
        >> q.mAlpha >> q.mDomainLength >> q.mNumBins >> q.mLambda
        >> q.mUseCachedRng >> q.mU1 >> q.mU2;
    return ar;
}