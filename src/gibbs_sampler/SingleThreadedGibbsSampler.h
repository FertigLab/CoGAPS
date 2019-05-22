#ifndef __COGAPS_SINGLE_THREADED_GIBBS_SAMPLER_H__
#define __COGAPS_SINGLE_THREADED_GIBBS_SAMPLER_H__

#include "AlphaParameters.h"
#include "../atomic/AtomicDomain.h"
#include "../atomic/ProposalQueue.h"
#include "../data_structures/Matrix.h"
#include "../math/Math.h"
#include "../math/VectorMath.h"
#include "../math/MatrixMath.h"
#include "../math/Random.h"
#include "../GapsParameters.h"

#include <vector>
#include <limits>
#include <cmath>
#include <stdint.h>

//////////////////////////// SingleThreadedGibbsSampler Interface ////////////////////////////

class GapsStatistics;

template <class StoragePolicy>
class SingleThreadedGibbsSampler;

template <class StoragePolicy>
Archive& operator<<(Archive &ar, const SingleThreadedGibbsSampler<StoragePolicy> &sampler);

template <class StoragePolicy>
Archive& operator>>(Archive &ar, SingleThreadedGibbsSampler<StoragePolicy> &sampler);

template <class StoragePolicy>
class SingleThreadedGibbsSampler : public StoragePolicy
{
public:

    friend class GapsStatistics;

    template <class DataType>
    SingleThreadedGibbsSampler(const DataType &data, bool transpose, bool subsetRows,
        float alpha, float maxGibbsMass, const GapsParameters &params,
        GapsRandomState *randState);

    void setAnnealingTemp(float temp);
    void setMatrix(const Matrix &mat);

    unsigned nAtoms() const;
    float dataSparsity() const;
    float getAverageQueueLength() const { return 0.f; }

    void update(unsigned nSteps, unsigned nThreads);

    friend Archive& operator<< <StoragePolicy> (Archive &ar, const SingleThreadedGibbsSampler &s);
    friend Archive& operator>> <StoragePolicy> (Archive &ar, SingleThreadedGibbsSampler &s);

#ifdef GAPS_DEBUG
    bool internallyConsistent() const;
#endif    

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    AtomicDomain mDomain; // data structure providing access to atoms

    GapsRandomState *mRandState;
    mutable GapsRng mRng;

    float mAlpha;
    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumPatterns;
    uint64_t mNumBins;
    uint64_t mBinLength;

    double mDomainLength; // length of entire atomic domain

    unsigned mNumCols;

    char getUpdateType();
    void birth();
    void death();
    void move();
    void exchange();

    bool canUseGibbs(unsigned col) const;
    bool canUseGibbs(unsigned c1, unsigned c2) const;
};

//////////////////// SingleThreadedGibbsSampler - templated functions ////////////////////////

template <class StoragePolicy>
template <class DataType>
SingleThreadedGibbsSampler<StoragePolicy>::SingleThreadedGibbsSampler(const DataType &data, bool transpose,
bool subsetRows, float alpha, float maxGibbsMass, const GapsParameters &params,
GapsRandomState *randState)
    :
StoragePolicy(data, transpose, subsetRows, params),
mDomain(StoragePolicy::mMatrix.nRow() * params.nPatterns),
mRandState(randState),
mRng(randState),
mAlpha(alpha),
mLambda(0.f),
mMaxGibbsMass(0.f),
mAnnealingTemp(1.f),
mNumPatterns(params.nPatterns),
mNumBins(StoragePolicy::mMatrix.nRow() * params.nPatterns),
mBinLength(std::numeric_limits<uint64_t>::max() / (StoragePolicy::mMatrix.nRow() * params.nPatterns)),
mDomainLength(static_cast<double>(mBinLength * static_cast<uint64_t>(StoragePolicy::mMatrix.nRow() * params.nPatterns))),
mNumCols(StoragePolicy::mMatrix.nCol())
{
    float meanD = params.singleCell ? gaps::nonZeroMean(StoragePolicy::mDMatrix) :
       gaps::mean(StoragePolicy::mDMatrix);

    mLambda = mAlpha * std::sqrt(mNumPatterns / meanD);
    mMaxGibbsMass = maxGibbsMass / mLambda;

    static bool warningGiven = false; // prevent duplicate warnings
    if (!warningGiven && gaps::max(StoragePolicy::mDMatrix) > 50.f)
    {
        warningGiven = true;
        gaps_printf("\nWarning: Large values detected in data, "
            "data needs to be log-transformed\n");
    }
}

template <class StoragePolicy>
void SingleThreadedGibbsSampler<StoragePolicy>::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

template <class StoragePolicy>
void SingleThreadedGibbsSampler<StoragePolicy>::setMatrix(const Matrix &mat)
{
    StoragePolicy::mMatrix = mat;
}

template <class StoragePolicy>
unsigned SingleThreadedGibbsSampler<StoragePolicy>::nAtoms() const
{
    return mDomain.size();
}

template <class StoragePolicy>
float SingleThreadedGibbsSampler<StoragePolicy>::dataSparsity() const
{
    return gaps::sparsity(StoragePolicy::mDMatrix);
}

template <class StoragePolicy>
char SingleThreadedGibbsSampler<StoragePolicy>::getUpdateType()
{
    if (mDomain.size() < 2)
    {
        return 'B'; // always birth when no atoms exist
    }

    double nAtoms = static_cast<double>(mDomain.size());
    double numer = nAtoms * mDomainLength;
    float deathProb = numer / (numer + mAlpha * mNumBins * (mDomainLength - nAtoms));

    float u1 = mRng.uniform();
    if (u1 < 0.5f)
    {
        return mRng.uniform() < deathProb ? 'D' : 'B';
    }
    return u1 < 0.75f ? 'M' : 'E';
}

template <class StoragePolicy>
void SingleThreadedGibbsSampler<StoragePolicy>::update(unsigned nSteps, unsigned nThreads)
{
    for (unsigned i = 0; i < nSteps; ++i)
    {
        switch (getUpdateType())
        {
            case 'B': birth();    break;
            case 'D': death();    break;
            case 'M': move();     break;
            case 'E': exchange(); break;
        }
    }

    //GAPS_ASSERT(internallyConsistent());
    GAPS_ASSERT(mDomain.isSorted());
}

inline float getDeltaLL_2(AlphaParameters alpha, float mass)
{
    return mass * (alpha.s_mu - alpha.s * mass / 2.f);
}

inline OptionalFloat gibbsMass_2(AlphaParameters alpha, float a, float b,
float lambda, GapsRng *rng)
{
    if (alpha.s > gaps::epsilon)
    {
        float mean = (alpha.s_mu - lambda) / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        return rng->truncNormal(a, b, mean, sd);
    }
    return OptionalFloat();
}

// add an atom at a random position, calculate mass either with an
// exponential distribution or with the gibbs mass distribution
template <class StoragePolicy>
void SingleThreadedGibbsSampler<StoragePolicy>::birth()
{
    uint64_t pos = mDomain.randomFreePosition(&mRng);
    unsigned row = (pos / mBinLength) / mNumCols;
    unsigned col = (pos / mBinLength) % mNumCols;

    // try to get mass using gibbs, resort to exponential if needed
    OptionalFloat mass = canUseGibbs(col)
        ? gibbsMass_2(StoragePolicy::alphaParameters(row, col) *
            mAnnealingTemp, 0.f, mMaxGibbsMass, mLambda, &mRng)
        : mRng.exponential(mLambda);

    // accept mass as long as gibbs succeded or it's non-zero
    if (mass.hasValue() && mass.value() >= gaps::epsilon)
    {
        mDomain.insert(pos, mass.value());
        StoragePolicy::changeMatrix(row, col, mass.value());
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass distribution
template <class StoragePolicy>
void SingleThreadedGibbsSampler<StoragePolicy>::death()
{
    Atom *atom = mDomain.randomAtom(&mRng);
    unsigned row = (atom->pos / mBinLength) / mNumCols;
    unsigned col = (atom->pos / mBinLength) % mNumCols;

    // calculate alpha parameters assuming atom dies
    AlphaParameters alpha = StoragePolicy::alphaParametersWithChange(row,
        col, -atom->mass);

    float rebirthMass = atom->mass;
    if (canUseGibbs(col))
    {
        OptionalFloat gMass = gibbsMass_2(alpha * mAnnealingTemp, 0.f,
            mMaxGibbsMass, mLambda, &mRng);
        if (gMass.hasValue())
        {
            rebirthMass = gMass.value();
        }
    }

    // accept/reject rebirth
    float deltaLL = getDeltaLL_2(alpha, rebirthMass) * mAnnealingTemp;
    if (std::log(mRng.uniform()) < deltaLL) // accept rebirth
    {
        if (rebirthMass != atom->mass)
        {
            StoragePolicy::safelyChangeMatrix(row, col,
                rebirthMass - atom->mass);
            atom->mass = rebirthMass;
        }
    }
    else // reject rebirth
    {
        StoragePolicy::safelyChangeMatrix(row, col, -atom->mass);
        mDomain.erase(atom->pos);
    }
}

// move mass from src to dest in the atomic domain
template <class StoragePolicy>
void SingleThreadedGibbsSampler<StoragePolicy>::move()
{
    AtomNeighborhood hood = mDomain.randomAtomWithNeighbors(&mRng);
    uint64_t lbound = hood.hasLeft() ? hood.left->pos : 0;
    uint64_t rbound = hood.hasRight() ? hood.right->pos : static_cast<uint64_t>(mDomainLength);

    uint64_t pos = mRng.uniform64(lbound + 1, rbound - 1);
    Atom *atom = hood.center;
    unsigned r1 = (atom->pos / mBinLength) / mNumCols;
    unsigned c1 = (atom->pos / mBinLength) % mNumCols;
    unsigned r2 = (pos / mBinLength) / mNumCols;
    unsigned c2 = (pos / mBinLength) % mNumCols;

    if (r1 == r2 && c1 == c2)
    {
        atom->pos = pos;
    }
    else
    {
        AlphaParameters alpha = StoragePolicy::alphaParameters(r1, c1, r2, c2);
        if (std::log(mRng.uniform()) < getDeltaLL_2(alpha, -atom->mass) * mAnnealingTemp)
        {
            atom->pos = pos;
            StoragePolicy::safelyChangeMatrix(r1, c1, -atom->mass);
            StoragePolicy::changeMatrix(r2, c2, atom->mass);
        }
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small
template <class StoragePolicy>
void SingleThreadedGibbsSampler<StoragePolicy>::exchange()
{
    AtomNeighborhood hood = mDomain.randomAtomWithRightNeighbor(&mRng);
    Atom *atom1 = hood.center;
    Atom *atom2 = hood.hasRight() ? hood.right : mDomain.front();
    unsigned r1 = (atom1->pos / mBinLength) / mNumCols;
    unsigned c1 = (atom1->pos / mBinLength) % mNumCols;
    unsigned r2 = (atom2->pos / mBinLength) / mNumCols;
    unsigned c2 = (atom2->pos / mBinLength) % mNumCols;

    if (r1 == r2 && c1 == c2)
    {
        float newMass = mRng.truncGammaUpper(atom1->mass + atom2->mass, 1.f / mLambda);
        float delta = (atom1->mass > atom2->mass) ? newMass - atom1->mass : atom2->mass - newMass;
        if (atom1->mass + delta > gaps::epsilon && atom2->mass - delta > gaps::epsilon)
        {
            atom1->mass += delta;
            atom2->mass -= delta;
        }
        return;
    }

    AlphaParameters alpha = StoragePolicy::alphaParameters(r1, c1, r2, c2);
    if (canUseGibbs(c1, c2))
    {
        OptionalFloat gMass = gibbsMass_2(alpha * mAnnealingTemp,
            -atom1->mass, atom2->mass, 0.f, &mRng);
        if (gMass.hasValue()
        && atom1->mass + gMass.value() > gaps::epsilon
        && atom2->mass - gMass.value() > gaps::epsilon)
        {
            atom1->mass += gMass.value();
            atom2->mass -= gMass.value();
            StoragePolicy::safelyChangeMatrix(r1, c1, gMass.value());
            StoragePolicy::safelyChangeMatrix(r2, c2, -gMass.value());
        }
        return;
    }
}

template <class StoragePolicy>
bool SingleThreadedGibbsSampler<StoragePolicy>::canUseGibbs(unsigned col) const
{
    return !gaps::isVectorZero(StoragePolicy::mOtherMatrix->getCol(col));
}

template <class StoragePolicy>
bool SingleThreadedGibbsSampler<StoragePolicy>::canUseGibbs(unsigned c1, unsigned c2) const
{
    return canUseGibbs(c1) || canUseGibbs(c2);
}

template <class StoragePolicy>
Archive& operator<<(Archive &ar, const SingleThreadedGibbsSampler<StoragePolicy> &s)
{
    operator<<(ar, static_cast<const StoragePolicy&>(s)) << s.mDomain
        << s.mAlpha << s.mLambda << s.mMaxGibbsMass
        << s.mAnnealingTemp << s.mNumPatterns << s.mNumBins << s.mBinLength;
    return ar;
}

template <class StoragePolicy>
Archive& operator>>(Archive &ar, SingleThreadedGibbsSampler<StoragePolicy> &s)
{
    operator>>(ar, static_cast<StoragePolicy&>(s)) >> s.mDomain
        >> s.mAlpha >> s.mLambda >> s.mMaxGibbsMass >> s.mAnnealingTemp
        >> s.mNumPatterns >> s.mNumBins >> s.mBinLength;
    return ar;
}

#ifdef GAPS_DEBUG
template <class StoragePolicy>
bool SingleThreadedGibbsSampler<StoragePolicy>::internallyConsistent() const
{
    return true;
}
#endif

#endif // __COGAPS_GIBBS_SAMPLER_H__