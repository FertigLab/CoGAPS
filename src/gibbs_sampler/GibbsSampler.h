#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#include "../atomic/AtomicDomain.h"
#include "../atomic/ProposalQueue.h"
#include "../data_structures/Matrix.h"
#include "../math/Math.h"
#include "../math/VectorMath.h"
#include "../math/MatrixMath.h"
#include "../math/Random.h"
#include "../GapsParameters.h"

#include <vector>

class GapsStatistics;

struct AlphaParameters
{
    float s;
    float s_mu;
    
    AlphaParameters(float t_s, float t_smu) : s(t_s), s_mu(t_smu) {}

    AlphaParameters operator+(const AlphaParameters &other) const
        { return AlphaParameters(s + other.s, s_mu - other.s_mu); } // not a typo

    AlphaParameters operator*(float v) const
        { return AlphaParameters(s * v, s_mu * v); }

    void operator*=(float v)
        { s *= v; s_mu *= v; }
};

//////////////////////// Common GibbsSampler Interface /////////////////////////

template <class Derived, class DataMatrix, class FactorMatrix>
class GibbsSampler
{
public:

    friend class GapsStatistics;

    void setAnnealingTemp(float temp);
    void setMatrix(const Matrix &mat);

    unsigned nAtoms() const;

    void update(unsigned nSteps, unsigned nThreads);

#ifndef GAPS_INTERNAL_TESTS
protected:
#endif

    DataMatrix mDMatrix; // samples by genes for A, genes by samples for P
    FactorMatrix mMatrix; // genes by patterns for A, samples by patterns for P
    const FactorMatrix *mOtherMatrix; // pointer to P if this is A, and vice versa

    AtomicDomain mDomain; // data structure providing access to atoms
    ProposalQueue mQueue; // creates queue of proposals that get evaluated by sampler

    float mAlpha;
    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumPatterns;
    uint64_t mNumBins;
    uint64_t mBinLength;

    // can't be constructed outside of derived classes
    template <class DataType>
    GibbsSampler(const DataType &data, bool transpose, bool subsetRows,
        float alpha, float maxGibbsMass, const GapsParameters &params,
        GapsRandomState *randState);

    Derived* impl();

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

template <class Derived, class DataMatrix, class FactorMatrix>
template <class DataType>
GibbsSampler<Derived, DataMatrix, FactorMatrix>::GibbsSampler(const DataType &data,
bool transpose, bool subsetRows, float alpha, float maxGibbsMass,
const GapsParameters &params, GapsRandomState *randState)
    :
mDMatrix(data, transpose, subsetRows, params.dataIndicesSubset),
mMatrix(mDMatrix.nCol(), params.nPatterns),
mOtherMatrix(NULL),
mDomain(mMatrix.nRow() * params.nPatterns),
mQueue(mMatrix.nRow(), params.nPatterns, randState),
mAlpha(alpha),
mLambda(0.f),
mMaxGibbsMass(0.f),
mAnnealingTemp(1.f),
mNumPatterns(params.nPatterns),
mNumBins(mMatrix.nRow() * params.nPatterns),
mBinLength(std::numeric_limits<uint64_t>::max() / mNumBins)
{
    float meanD = params.singleCell ? gaps::nonZeroMean(mDMatrix) :
       gaps::mean(mDMatrix);

    mLambda = mAlpha * std::sqrt(mNumPatterns / meanD);
    mMaxGibbsMass = maxGibbsMass / mLambda;
    mQueue.setAlpha(mAlpha);
    mQueue.setLambda(mLambda);
}

template <class Derived, class DataMatrix, class FactorMatrix>
Derived* GibbsSampler<Derived, DataMatrix, FactorMatrix>::impl()
{
    return static_cast<Derived*>(this);
}

template <class Derived, class DataMatrix, class FactorMatrix>
void GibbsSampler<Derived, DataMatrix, FactorMatrix>::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

template <class Derived, class DataMatrix, class FactorMatrix>
void GibbsSampler<Derived, DataMatrix, FactorMatrix>::setMatrix(const Matrix &mat)
{
    mMatrix = mat;
}

template <class Derived, class DataMatrix, class FactorMatrix>
unsigned GibbsSampler<Derived, DataMatrix, FactorMatrix>::nAtoms() const
{
    return mDomain.size();
}

template <class Derived, class DataMatrix, class FactorMatrix>
void GibbsSampler<Derived, DataMatrix, FactorMatrix>::update(unsigned nSteps, unsigned nThreads)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        // create the largest queue possible, without hitting any conflicts
        mQueue.populate(mDomain, nSteps - n);
        n += mQueue.size();
        
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

    GAPS_ASSERT(impl()->internallyConsistent());
    GAPS_ASSERT(mDomain.isSorted());
}

static float getDeltaLL(AlphaParameters alpha, float mass)
{
    return mass * (alpha.s_mu - alpha.s * mass / 2.f);
}

static OptionalFloat gibbsMass(AlphaParameters alpha, float a, float b,
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
template <class Derived, class DataMatrix, class FactorMatrix>
void GibbsSampler<Derived, DataMatrix, FactorMatrix>::birth(const AtomicProposal &prop)
{
    // try to get mass using gibbs, resort to exponential if needed
    OptionalFloat mass = canUseGibbs(prop.c1)
        ? gibbsMass(impl()->alphaParameters(prop.r1, prop.c1) * mAnnealingTemp,
            0.f, mMaxGibbsMass, mLambda, &(prop.rng))
        : prop.rng.exponential(mLambda);

    // accept mass as long as gibbs succeded or it's non-zero
    if (mass.hasValue() && mass.value() >= gaps::epsilon)
    {
        mQueue.acceptBirth();
        prop.atom1->mass = mass.value();
        impl()->changeMatrix(prop.r1, prop.c1, mass.value());
    }
    else
    {
        mQueue.rejectBirth();
        mDomain.erase(prop.atom1->pos);
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass distribution
template <class Derived, class DataMatrix, class FactorMatrix>
void GibbsSampler<Derived, DataMatrix, FactorMatrix>::death(const AtomicProposal &prop)
{
    // calculate alpha parameters assuming atom dies
    AlphaParameters alpha = impl()->alphaParametersWithChange(prop.r1, prop.c1,
        -prop.atom1->mass);

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
            impl()->safelyChangeMatrix(prop.r1, prop.c1, rebirthMass - prop.atom1->mass);
        }
        prop.atom1->mass = rebirthMass;
    }
    else // reject rebirth
    {
        mQueue.acceptDeath();
        impl()->safelyChangeMatrix(prop.r1, prop.c1, -prop.atom1->mass);
        mDomain.erase(prop.atom1->pos);
    }
}

// move mass from src to dest in the atomic domain
template <class Derived, class DataMatrix, class FactorMatrix>
void GibbsSampler<Derived, DataMatrix, FactorMatrix>::move(const AtomicProposal &prop)
{
    GAPS_ASSERT(prop.r1 != prop.r2 || prop.c1 != prop.c2);

    AlphaParameters alpha = impl()->alphaParameters(prop.r1, prop.c1, prop.r2, prop.c2);
    if (std::log(prop.rng.uniform()) < getDeltaLL(alpha, -prop.atom1->mass) * mAnnealingTemp)
    {
        prop.atom1->pos = prop.pos;
        impl()->safelyChangeMatrix(prop.r1, prop.c1, -prop.atom1->mass);
        impl()->changeMatrix(prop.r2, prop.c2, prop.atom1->mass);
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small
template <class Derived, class DataMatrix, class FactorMatrix>
void GibbsSampler<Derived, DataMatrix, FactorMatrix>::exchange(const AtomicProposal &prop)
{
    GAPS_ASSERT(prop.r1 != prop.r2 || prop.c1 != prop.c2);

    // attempt gibbs distribution exchange
    AlphaParameters alpha = impl()->alphaParameters(prop.r1, prop.c1, prop.r2, prop.c2);
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

template <class Derived, class DataMatrix, class FactorMatrix>
void GibbsSampler<Derived, DataMatrix, FactorMatrix>::exchangeUsingMetropolisHastings(const AtomicProposal &prop,
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
template <class Derived, class DataMatrix, class FactorMatrix>
void GibbsSampler<Derived, DataMatrix, FactorMatrix>::acceptExchange(const AtomicProposal &prop,
float delta)
{
    if (prop.atom1->mass + delta > gaps::epsilon && prop.atom2->mass - delta > gaps::epsilon)
    {
        prop.atom1->mass += delta;
        prop.atom2->mass -= delta;

        impl()->safelyChangeMatrix(prop.r1, prop.c1, delta);
        impl()->safelyChangeMatrix(prop.r2, prop.c2, -delta);
    }
}

template <class Derived, class DataMatrix, class FactorMatrix>
bool GibbsSampler<Derived, DataMatrix, FactorMatrix>::canUseGibbs(unsigned col) const
{
    return !gaps::isVectorZero(mMatrix.getCol(col));
}

template <class Derived, class DataMatrix, class FactorMatrix>
bool GibbsSampler<Derived, DataMatrix, FactorMatrix>::canUseGibbs(unsigned c1, unsigned c2) const
{
    return canUseGibbs(c1) || canUseGibbs(c2);
}

#endif // __COGAPS_GIBBS_SAMPLER_H__