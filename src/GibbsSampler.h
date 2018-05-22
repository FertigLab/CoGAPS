#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#include "GapsAssert.h"
#include "Archive.h"
#include "data_structures/Matrix.h"
#include "math/Random.h"
#include "math/Algorithms.h"
#include "ProposalQueue.h"
#include "AtomicDomain.h"

#include <Rcpp.h>
#include <algorithm>
#include <omp.h>

// forward declarations needed for friend classes
class AmplitudeGibbsSampler;
class PatternGibbsSampler;
class GapsStatistics;

/************************** GIBBS SAMPLER INTERFACE **************************/

template <class T, class MatA, class MatB>
class GibbsSampler
{
private:

    friend T; // prevent incorrect inheritance - only T can construct
    friend GapsStatistics;

    GibbsSampler(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
        unsigned nrow, unsigned ncol, float alpha);

protected:

    MatA mMatrix;
    MatB* mOtherMatrix;
    MatB mDMatrix;
    MatB mSMatrix;
    MatB mAPMatrix;

    ProposalQueue mQueue;
    AtomicDomain mDomain;

    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumRows;
    unsigned mNumCols;
    uint64_t mBinSize;

    float mAvgQueue;
    float mNumQueues;

    T* impl();

    void processProposal(const AtomicProposal &prop);

    void birth(uint64_t pos, unsigned row, unsigned col);
    void death(uint64_t pos, float mass, unsigned row, unsigned col);
    void move(uint64_t src, float mass, uint64_t dest, unsigned r1, unsigned c1,
        unsigned r2, unsigned c2);
    void exchange(uint64_t p1, float mass1, uint64_t p2, float mass2,
        unsigned r1, unsigned c1, unsigned r2, unsigned c2);

    float gibbsMass(unsigned row, unsigned col, float mass);
    float gibbsMass(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);

    void addMass(uint64_t pos, float mass, unsigned row, unsigned col);
    void removeMass(uint64_t pos, float mass, unsigned row, unsigned col);
    float updateAtomMass(uint64_t pos, float mass, float delta);

    void acceptExchange(uint64_t p1, float m1, float d1, uint64_t p2, float m2,
        float d2, unsigned r1, unsigned c1, unsigned r2, unsigned c2);

public:

    void update(unsigned nSteps, unsigned nCores);
    void setAnnealingTemp(float temp);
    float getAvgQueue() const { return mAvgQueue; }
    
    float chi2() const;
    uint64_t nAtoms() const;
};

class AmplitudeGibbsSampler : public GibbsSampler<AmplitudeGibbsSampler, ColMatrix, RowMatrix>
{
private:

    friend GibbsSampler;
    friend PatternGibbsSampler;

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;
    bool canUseGibbs(unsigned row, unsigned col) const;
    bool canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const;
    void updateAPMatrix(unsigned row, unsigned col, float delta);

public:

    AmplitudeGibbsSampler(const Rcpp::NumericMatrix &D,
        const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha=0.f,
        float maxGibbsMass=0.f);

    void sync(PatternGibbsSampler &sampler);

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2,
        unsigned c2);

    float computeDeltaLL(unsigned row, unsigned col, float mass);
    float computeDeltaLL(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);
};

class PatternGibbsSampler : public GibbsSampler<PatternGibbsSampler, RowMatrix, ColMatrix>
{
private:

    friend GibbsSampler;
    friend AmplitudeGibbsSampler;

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;
    bool canUseGibbs(unsigned row, unsigned col) const;
    bool canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const;
    void updateAPMatrix(unsigned row, unsigned col, float delta);

public:

    PatternGibbsSampler(const Rcpp::NumericMatrix &D,
        const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha=0.f,
        float maxGibbsMass=0.f);

    void sync(AmplitudeGibbsSampler &sampler);

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2,
        unsigned c2);

    float computeDeltaLL(unsigned row, unsigned col, float mass);
    float computeDeltaLL(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);
};

/******************* IMPLEMENTATION OF TEMPLATED FUNCTIONS *******************/

template <class T, class MatA, class MatB>
GibbsSampler<T, MatA, MatB>::GibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nrow, unsigned ncol, float alpha)
    :
mMatrix(nrow, ncol), mOtherMatrix(NULL), mDMatrix(D), mSMatrix(S),
mAPMatrix(D.nrow(), D.ncol()), mQueue(nrow * ncol, alpha),
mAnnealingTemp(0.f), mNumRows(nrow), mNumCols(ncol),
mAvgQueue(0.f), mNumQueues(0.f)
{
    mBinSize = std::numeric_limits<uint64_t>::max()
        / static_cast<uint64_t>(mNumRows * mNumCols);
    uint64_t remain = std::numeric_limits<uint64_t>::max()
        % static_cast<uint64_t>(mNumRows * mNumCols);
    mQueue.setDomainSize(std::numeric_limits<uint64_t>::max() - remain);
    mDomain.setDomainSize(std::numeric_limits<uint64_t>::max() - remain);
}

template <class T, class MatA, class MatB>
T* GibbsSampler<T, MatA, MatB>::impl()
{
    return static_cast<T*>(this);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::update(unsigned nSteps, unsigned nCores)
{
    static unsigned count = 0;
    unsigned n = 0;
    while (n < nSteps)
    {
        //GAPS_ASSERT(nSteps - (mQueue.size() + n) >= 0);
        mQueue.populate(mDomain, nSteps - n);
        //mQueue.populate(mDomain, nSteps - (mQueue.size() + n));
        
        //GAPS_ASSERT(mQueue.size() > 0);
        mNumQueues += 1.f;
        mAvgQueue = mQueue.size() / mNumQueues + mAvgQueue * (mNumQueues - 1.f) / mNumQueues;
        n += mQueue.size();
        mDomain.resetCache(mQueue.size());
        //Rprintf("round: %d\n", count++);

	//unsigned nThreads = std::min(mQueue.size(), nCores);
        #pragma omp parallel for num_threads(nCores)
        for (unsigned i = 0; i < mQueue.size(); ++i)
        {
            processProposal(mQueue[i]);
        }
        mDomain.flushCache();
        mQueue.clear(1);
        //GAPS_ASSERT(n <= nSteps);
    }
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::processProposal(const AtomicProposal &prop)
{
    unsigned r1 = impl()->getRow(prop.pos1);
    unsigned c1 = impl()->getCol(prop.pos1);
    unsigned r2 = 0, c2 = 0;
    //Rprintf("type: %c pos1: %lu\n", prop.type, prop.pos1);
    switch (prop.type)
    {
        case 'B':
            birth(prop.pos1, r1, c1);
            break;
        case 'D':
            death(prop.pos1, prop.mass1, r1, c1);
            break;
        case 'M':
            r2 = impl()->getRow(prop.pos2);
            c2 = impl()->getCol(prop.pos2);
            move(prop.pos1, prop.mass1, prop.pos2, r1, c1, r2, c2);
            break;
        case 'E':
            //Rprintf("pos2: %lu\n", prop.pos2);
            r2 = impl()->getRow(prop.pos2);
            c2 = impl()->getCol(prop.pos2);
            exchange(prop.pos1, prop.mass1, prop.pos2, prop.mass2, r1, c1, r2, c2);
            break;
    }
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::addMass(uint64_t pos, float mass, unsigned row, unsigned col)
{
    mDomain.cacheInsert(pos, mass);
    mMatrix(row, col) += mass;
    impl()->updateAPMatrix(row, col, mass);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::removeMass(uint64_t pos, float mass, unsigned row, unsigned col)
{
    mDomain.cacheErase(pos);
    mMatrix(row, col) += -mass;
    impl()->updateAPMatrix(row, col, -mass);
}

// add an atom at pos, calculate mass either with an exponential distribution
// or with the gibbs mass calculation
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::birth(uint64_t pos, unsigned row,
unsigned col)
{
    float mass = impl()->canUseGibbs(row, col) ? gibbsMass(row, col, mass)
        : gaps::random::exponential(mLambda);
    if (mass >= gaps::algo::epsilon)
    {
        mDomain.updateMass(pos, mass);
        mMatrix(row, col) += mass;
        impl()->updateAPMatrix(row, col, mass);
    }
    else
    {
        mDomain.cacheErase(pos);
        mQueue.rejectBirth();
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass calculation
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::death(uint64_t pos, float mass, unsigned row,
unsigned col)
{
    //GAPS_ASSERT(mass > 0.f);

    //removeMass(pos, mass, row, col);
    mMatrix(row, col) += -mass;
    impl()->updateAPMatrix(row, col, -mass);

    float newMass = impl()->canUseGibbs(row, col) ? gibbsMass(row, col, -mass) : 0.f;
    mass = newMass < gaps::algo::epsilon ? mass : newMass;
    float deltaLL = impl()->computeDeltaLL(row, col, mass);
    if (deltaLL * mAnnealingTemp >= std::log(gaps::random::uniform()))
    {
        mDomain.updateMass(pos, mass);
        mMatrix(row, col) += mass;
        impl()->updateAPMatrix(row, col, mass);
        mQueue.rejectDeath();
    }
    else
    {
        mDomain.cacheErase(pos);
        mQueue.acceptDeath();
    }
}

// move mass from src to dest in the atomic domain
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::move(uint64_t src, float mass, uint64_t dest,
unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    //GAPS_ASSERT(mass > 0.f);
    if (r1 != r2 || c1 != c2) // automatically reject if change in same bin
    {
        float deltaLL = impl()->computeDeltaLL(r1, c1, -mass, r2, c2, mass);
        if (deltaLL * mAnnealingTemp > std::log(gaps::random::uniform()))
        {
            removeMass(src, mass, r1, c1);
            addMass(dest, mass, r2, c2);
        }
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small after
// the exchange
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::exchange(uint64_t p1, float m1, uint64_t p2,
float m2, unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    //GAPS_ASSERT(m1 > 0.f);
    //GAPS_ASSERT(m2 > 0.f);
    float pUpper = gaps::random::p_gamma(m1 + m2, 2.f, 1.f / mLambda);
    float newMass = gaps::random::inverseGammaSample(0.f, pUpper, 2.f, 1.f / mLambda);
    if (r1 != r2 || c1 != c2) // automatically reject if change in same bin
    {
        if ((m1 > m2 && newMass - m1 < 0) || (m1 < m2 && m2 - newMass < 0))
        {
            std::swap(r1, r2);
            std::swap(c1, c2);
            std::swap(p1, p2);
            std::swap(m1, m2);
        }

        if (impl()->canUseGibbs(r1, c1, r2, c2))
        {
            float gDelta = gibbsMass(r1, c1, m1, r2, c2, m2);
            if (gDelta > -m1 - 0.5f) // janky, should be pair<bool, float>
            {
                acceptExchange(p1, m1, gDelta, p2, m2, -gDelta, r1, c1, r2, c2);
                return;
            }
        }

        // use metropolis hastings otherwise
        float delta = m1 > m2 ? newMass - m1 : m2 - newMass; // change larger mass
        float pOldMass = m1 + delta > m2 - delta ? m1 : m2;
        float pNew = gaps::random::d_gamma(newMass, 2.f, 1.f / mLambda);
        float pOld = gaps::random::d_gamma(pOldMass, 2.f, 1.f / mLambda);

        if (pOld == 0.f && pNew != 0.f) // special case
        {
            acceptExchange(p1, m1, delta, p2, m2, -delta, r1, c1, r2, c2);
            return;
        }
        else // normal case
        {
            float deltaLL = impl()->computeDeltaLL(r1, c1, delta, r2, c2, -delta);
            float priorLL = (pOld == 0.f) ? 1.f : pOld / pNew;
            float u = std::log(gaps::random::uniform() * priorLL);
            if (u < deltaLL * mAnnealingTemp)
            {
                acceptExchange(p1, m1, delta, p2, m2, -delta, r1, c1, r2, c2);
                return;
            }
        }
    }
    mQueue.rejectDeath();
}

// helper function for acceptExchange, used to conditionally update the mass
// at a single position, deleting the atom if the new mass is too small
template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::updateAtomMass(uint64_t pos, float mass,
float delta)
{
    if (mass + delta < gaps::algo::epsilon)
    {
        mDomain.cacheErase(pos);
        mQueue.acceptDeath();
        return false;
    }
    else
    {
        mDomain.updateMass(pos, mass + delta);
        return true;
    }
}

// helper function for exchange step, updates the atomic domain, matrix, and
// A*P matrix
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::acceptExchange(uint64_t p1, float m1,
float d1, uint64_t p2, float m2, float d2, unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    //Rprintf("p1: %lu p2: %lu\n", p1, p2);
    bool b1 = updateAtomMass(p1, m1, d1);
    bool b2 = updateAtomMass(p2, m2, d2);
    //GAPS_ASSERT(b1 || b2);
    
    if (!b1) d1 = -m1;
    if (!b2) d2 = -m2;

    if (b1 && b2)
    {
        mQueue.rejectDeath();
    }

    mMatrix(r1, c1) += d1;
    mMatrix(r2, c2) += d2;
    impl()->updateAPMatrix(r1, c1, d1);
    impl()->updateAPMatrix(r2, c2, d2);
}

template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::gibbsMass(unsigned row, unsigned col, float mass)
{        
    AlphaParameters alpha = impl()->alphaParameters(row, col);
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::algo::epsilon)
    {
        float mean = (alpha.su - mLambda) / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::random::p_norm(0.f, mean, sd);

        if (pLower < 1.f)
        {
            float m = gaps::random::inverseNormSample(pLower, 1.f, mean, sd);
            return std::max(std::min(m, mMaxGibbsMass), 0.f);
        }
    }
    return std::min(mass < 0.f ? std::abs(mass) : 0.f, mMaxGibbsMass);
}

template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::gibbsMass(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2)
{
    AlphaParameters alpha = impl()->alphaParameters(r1, c1, r2, c2);
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::algo::epsilon)
    {
        float mean = alpha.su / alpha.s; // TODO why not subtract lambda
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::random::p_norm(-m1, mean, sd);
        float pUpper = gaps::random::p_norm(m2, mean, sd);

        if (!(pLower >  0.95f || pUpper < 0.05f))
        {
            float delta = gaps::random::inverseNormSample(pLower, pUpper, mean, sd);
            return std::min(std::max(-m1, delta), m2); // conserve mass
        }
    }
    return -m1 - 1.f;
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}
  
template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::chi2() const
{   
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}
  
template <class T, class MatA, class MatB>
uint64_t GibbsSampler<T, MatA, MatB>::nAtoms() const
{   
    return mDomain.size();
}

#endif
