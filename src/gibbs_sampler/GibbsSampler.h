#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

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

//////////////////////////// GibbsSampler Interface ////////////////////////////

class GapsStatistics;

template <class StoragePolicy>
class GibbsSampler;

template <class StoragePolicy>
Archive& operator<<(Archive &ar, const GibbsSampler<StoragePolicy> &sampler);

template <class StoragePolicy>
Archive& operator>>(Archive &ar, GibbsSampler<StoragePolicy> &sampler);

template <class StoragePolicy>
class GibbsSampler : public StoragePolicy
{
public:

    friend class GapsStatistics;

    template <class DataType>
    GibbsSampler(const DataType &data, bool transpose, bool subsetRows,
        float alpha, float maxGibbsMass, const GapsParameters &params,
        GapsRandomState *randState);

    void setAnnealingTemp(float temp);
    void setMatrix(const Matrix &mat);

    unsigned nAtoms() const;

    void update(unsigned nSteps, unsigned nThreads);

    friend Archive& operator<< <StoragePolicy> (Archive &ar, const GibbsSampler &s);
    friend Archive& operator>> <StoragePolicy> (Archive &ar, GibbsSampler &s);

#ifdef GAPS_DEBUG
    bool internallyConsistent() const;
#endif    

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    AtomicDomain mDomain; // data structure providing access to atoms
    ProposalQueue mQueue; // creates queue of proposals that get evaluated by sampler

    float mAlpha;
    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumPatterns;
    uint64_t mNumBins;
    uint64_t mBinLength;

    void birth(const AtomicProposal &prop);
    void death(const AtomicProposal &prop);
    void move(const AtomicProposal &prop);
    void exchange(const AtomicProposal &prop);
    void exchangeUsingMetropolisHastings(const AtomicProposal &prop,
        AlphaParameters alpha);
    void acceptExchange(const AtomicProposal &prop, float delta);

    bool canUseGibbs(unsigned col) const;
    bool canUseGibbs(unsigned c1, unsigned c2) const;
};

//////////////////// GibbsSampler - templated functions ////////////////////////

template <class StoragePolicy>
template <class DataType>
GibbsSampler<StoragePolicy>::GibbsSampler(const DataType &data, bool transpose,
bool subsetRows, float alpha, float maxGibbsMass, const GapsParameters &params,
GapsRandomState *randState)
    :
StoragePolicy(data, transpose, subsetRows, params),
mDomain(StoragePolicy::mMatrix.nRow() * params.nPatterns),
mQueue(StoragePolicy::mMatrix.nRow(), params.nPatterns, randState),
mAlpha(alpha),
mLambda(0.f),
mMaxGibbsMass(0.f),
mAnnealingTemp(1.f),
mNumPatterns(params.nPatterns),
mNumBins(StoragePolicy::mMatrix.nRow() * params.nPatterns),
mBinLength(std::numeric_limits<uint64_t>::max() / (StoragePolicy::mMatrix.nRow() * params.nPatterns))
{
    float meanD = params.singleCell ? gaps::nonZeroMean(StoragePolicy::mDMatrix) :
       gaps::mean(StoragePolicy::mDMatrix);

    mLambda = mAlpha * std::sqrt(mNumPatterns / meanD);
    mMaxGibbsMass = maxGibbsMass / mLambda;
    mQueue.setAlpha(mAlpha);
    mQueue.setLambda(mLambda);

    static bool warningGiven = false; // prevent duplicate warnings
    if (!warningGiven && gaps::max(StoragePolicy::mDMatrix) > 50.f)
    {
        warningGiven = true;
        gaps_printf("\nWarning: Large values detected in data, "
            "data needs to be log-transformed\n");
    }
}

template <class StoragePolicy>
void GibbsSampler<StoragePolicy>::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

template <class StoragePolicy>
void GibbsSampler<StoragePolicy>::setMatrix(const Matrix &mat)
{
    StoragePolicy::mMatrix = mat;
}

template <class StoragePolicy>
unsigned GibbsSampler<StoragePolicy>::nAtoms() const
{
    return mDomain.size();
}

template <class StoragePolicy>
void GibbsSampler<StoragePolicy>::update(unsigned nSteps, unsigned nThreads)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        // create the largest queue possible, without hitting any conflicts
        mQueue.populate(mDomain, nSteps - n);
        n += mQueue.nProcessed();
        
        // process all proposed updates in parallel - the way the queue is 
        // populated ensures no race conditions will happen
        #pragma omp parallel for num_threads(nThreads)
        for (unsigned i = 0; i < mQueue.size(); ++i)
        {
            switch (mQueue[i].type)
            {
                case 'B': birth(mQueue[i]);    break;
                case 'D': death(mQueue[i]);    break;
                case 'M': move(mQueue[i]);     break;
                case 'E': exchange(mQueue[i]); break;
            }
        }
        mQueue.clear();
    }

    //GAPS_ASSERT(internallyConsistent());
    GAPS_ASSERT(mDomain.isSorted());
}

inline float getDeltaLL(AlphaParameters alpha, float mass)
{
    return mass * (alpha.s_mu - alpha.s * mass / 2.f);
}

inline OptionalFloat gibbsMass(AlphaParameters alpha, float a, float b,
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
void GibbsSampler<StoragePolicy>::birth(const AtomicProposal &prop)
{
    // try to get mass using gibbs, resort to exponential if needed
    OptionalFloat mass = canUseGibbs(prop.c1)
        ? gibbsMass(StoragePolicy::alphaParameters(prop.r1, prop.c1) *
            mAnnealingTemp, 0.f, mMaxGibbsMass, mLambda, &(prop.rng))
        : prop.rng.exponential(mLambda);

    // accept mass as long as gibbs succeded or it's non-zero
    if (mass.hasValue() && mass.value() >= gaps::epsilon)
    {
        mQueue.acceptBirth();
        prop.atom1->mass = mass.value();
        StoragePolicy::changeMatrix(prop.r1, prop.c1, mass.value());
    }
    else
    {
        mQueue.rejectBirth();
        mDomain.erase(prop.atom1->pos);
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass distribution
template <class StoragePolicy>
void GibbsSampler<StoragePolicy>::death(const AtomicProposal &prop)
{
    // calculate alpha parameters assuming atom dies
    AlphaParameters alpha = StoragePolicy::alphaParametersWithChange(prop.r1,
        prop.c1, -prop.atom1->mass);

    float rebirthMass = prop.atom1->mass;
    if (canUseGibbs(prop.c1))
    {
        OptionalFloat gMass = gibbsMass(alpha * mAnnealingTemp, 0.f,
            mMaxGibbsMass, mLambda, &(prop.rng));
        if (gMass.hasValue())
        {
            rebirthMass = gMass.value();
        }
    }

    // accept/reject rebirth
    float deltaLL = getDeltaLL(alpha, rebirthMass) * mAnnealingTemp;
    if (std::log(prop.rng.uniform()) < deltaLL) // accept rebirth
    {
        mQueue.rejectDeath();
        if (rebirthMass != prop.atom1->mass)
        {
            StoragePolicy::safelyChangeMatrix(prop.r1, prop.c1,
                rebirthMass - prop.atom1->mass);
            prop.atom1->mass = rebirthMass;
        }
    }
    else // reject rebirth
    {
        mQueue.acceptDeath();
        StoragePolicy::safelyChangeMatrix(prop.r1, prop.c1, -prop.atom1->mass);
        mDomain.erase(prop.atom1->pos);
    }
}

// move mass from src to dest in the atomic domain
template <class StoragePolicy>
void GibbsSampler<StoragePolicy>::move(const AtomicProposal &prop)
{
    GAPS_ASSERT(prop.r1 != prop.r2 || prop.c1 != prop.c2);

    AlphaParameters alpha = StoragePolicy::alphaParameters(prop.r1, prop.c1,
        prop.r2, prop.c2);
    if (std::log(prop.rng.uniform()) < getDeltaLL(alpha, -prop.atom1->mass) * mAnnealingTemp)
    {
        prop.atom1->pos = prop.pos;
        StoragePolicy::safelyChangeMatrix(prop.r1, prop.c1, -prop.atom1->mass);
        StoragePolicy::changeMatrix(prop.r2, prop.c2, prop.atom1->mass);
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small
template <class StoragePolicy>
void GibbsSampler<StoragePolicy>::exchange(const AtomicProposal &prop)
{
    GAPS_ASSERT(prop.r1 != prop.r2 || prop.c1 != prop.c2);

    // attempt gibbs distribution exchange
    AlphaParameters alpha = StoragePolicy::alphaParameters(prop.r1, prop.c1,
        prop.r2, prop.c2);
    if (canUseGibbs(prop.c1, prop.c2))
    {
        OptionalFloat gMass = gibbsMass(alpha * mAnnealingTemp,
            -prop.atom1->mass, prop.atom2->mass, 0.f, &(prop.rng));
        if (gMass.hasValue())
        {
            acceptExchange(prop, gMass.value());
            return;
        }
    }

    // resort to metropolis-hastings if gibbs fails
    exchangeUsingMetropolisHastings(prop, alpha);
}

template <class StoragePolicy>
void GibbsSampler<StoragePolicy>::exchangeUsingMetropolisHastings(const AtomicProposal &prop,
AlphaParameters alpha)
{
    // compute amount of mass to be exchanged
    float totalMass = prop.atom1->mass + prop.atom2->mass;
    float newMass = prop.rng.truncGammaUpper(totalMass, 1.f / mLambda);

    // compute amount to change atom1 by - always change larger mass to newMass
    float delta = (prop.atom1->mass > prop.atom2->mass)
        ? newMass - prop.atom1->mass
        : prop.atom2->mass - newMass;

    // choose mass for prior likelihood calculation
    float oldMass = (2.f * newMass > totalMass)
        ? gaps::max(prop.atom1->mass, prop.atom2->mass)
        : gaps::min(prop.atom1->mass, prop.atom2->mass);

    // calculate prior likelihood
    float pNew = gaps::d_gamma(newMass, 2.f, 1.f / mLambda);
    float pOld = gaps::d_gamma(oldMass, 2.f, 1.f / mLambda);
    float priorLL = (pNew == 0.f) ? 1.f : pOld / pNew;

    // accept/reject
    float deltaLL = getDeltaLL(alpha, delta) * mAnnealingTemp;
    if (priorLL == 0.f || std::log(prop.rng.uniform() * priorLL) < deltaLL)
    {
        acceptExchange(prop, delta);
        return;
    }
}

// helper function for exchange step
template <class StoragePolicy>
void GibbsSampler<StoragePolicy>::acceptExchange(const AtomicProposal &prop,
float delta)
{
    if (prop.atom1->mass + delta > gaps::epsilon
    && prop.atom2->mass - delta > gaps::epsilon)
    {
        float newMass1 = prop.atom1->mass + delta;
        float newMass2 = prop.atom2->mass - delta;
    
        StoragePolicy::safelyChangeMatrix(prop.r1, prop.c1,
            newMass1 - prop.atom1->mass);
        StoragePolicy::safelyChangeMatrix(prop.r2, prop.c2,
            newMass2 - prop.atom2->mass);

        prop.atom1->mass = newMass1;
        prop.atom2->mass = newMass2;
    }
}

template <class StoragePolicy>
bool GibbsSampler<StoragePolicy>::canUseGibbs(unsigned col) const
{
    return !gaps::isVectorZero(StoragePolicy::mOtherMatrix->getCol(col));
}

template <class StoragePolicy>
bool GibbsSampler<StoragePolicy>::canUseGibbs(unsigned c1, unsigned c2) const
{
    return canUseGibbs(c1) || canUseGibbs(c2);
}

template <class StoragePolicy>
Archive& operator<<(Archive &ar, const GibbsSampler<StoragePolicy> &s)
{
    operator<<(ar, static_cast<const StoragePolicy&>(s)) << s.mDomain
        << s.mQueue << s.mAlpha << s.mLambda << s.mMaxGibbsMass
        << s.mAnnealingTemp << s.mNumPatterns << s.mNumBins << s.mBinLength;
    return ar;
}

template <class StoragePolicy>
Archive& operator>>(Archive &ar, GibbsSampler<StoragePolicy> &s)
{
    operator>>(ar, static_cast<StoragePolicy&>(s)) >> s.mDomain >> s.mQueue
        >> s.mAlpha >> s.mLambda >> s.mMaxGibbsMass >> s.mAnnealingTemp
        >> s.mNumPatterns >> s.mNumBins >> s.mBinLength;
    return ar;
}

#ifdef GAPS_DEBUG
template <class StoragePolicy>
bool GibbsSampler<StoragePolicy>::internallyConsistent() const
{
    return true;
}
#endif

#endif // __COGAPS_GIBBS_SAMPLER_H__