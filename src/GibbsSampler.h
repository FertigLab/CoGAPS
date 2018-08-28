#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#include "AtomicDomain.h"
#include "ProposalQueue.h"
#include "data_structures/Matrix.h"
#include "math/Algorithms.h"
#include "math/Random.h"
#include "utils/Archive.h"
#include "utils/GapsAssert.h"

#include <algorithm>

// forward declarations needed for friend classes/functions
class AmplitudeGibbsSampler;
class PatternGibbsSampler;
class GapsStatistics;

template <class T, class MatA, class MatB>
class GibbsSampler;

template <class T, class MatA, class MatB>
inline Archive& operator<<(Archive &ar, GibbsSampler<T, MatA, MatB> &samp);

template <class T, class MatA, class MatB>
inline Archive& operator>>(Archive &ar, GibbsSampler<T, MatA, MatB> &samp);

/*************************** GIBBS SAMPLER INTERFACE **************************/

template <class T, class MatA, class MatB>
class GibbsSampler
{
private:

    friend class GapsStatistics;

protected:

    MatB mDMatrix;
    MatB mSMatrix;
    MatB mAPMatrix;

    ColMatrix mMatrix;
    const ColMatrix* mOtherMatrix;

    GapsRng mPropRng;
    AtomicDomain mDomain;

    float mAlpha;
    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumRows;
    unsigned mNumCols;

    uint64_t mNumBins;
    uint64_t mBinSize;
    uint64_t mDomainLength;

    float mAvgQueue;
    float mNumQueues;

    T* impl();

    void makeAndProcessProposal();
    float deathProb(uint64_t nAtoms) const;

    void birth();
    void death();
    void move();
    void exchange();

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;

    bool canUseGibbs(unsigned col) const;
    bool canUseGibbs(unsigned c1, unsigned c2) const;

    bool updateAtomMass(Atom *atom, float delta);
    void acceptExchange(AtomicProposal *prop, float d1, unsigned r1,
        unsigned c1, unsigned r2, unsigned c2);

    OptionalFloat gibbsMass(AlphaParameters alpha, GapsRng *rng);
    OptionalFloat gibbsMass(AlphaParameters alpha, float m1, float m2, GapsRng *rng);

public:

    template <class DataType>
    GibbsSampler(const DataType &data, bool transpose, unsigned nPatterns,
        bool amp, bool partitionRows, const std::vector<unsigned> &indices);

    template <class DataType>
    void setUncertainty(const DataType &unc, bool transposeData,
        bool partitionRows, const std::vector<unsigned> &indices);
    
    void setSparsity(float alpha, bool singleCell);
    void setMaxGibbsMass(float max);
    void setAnnealingTemp(float temp);

    void setMatrix(const Matrix &mat);

    void update(unsigned nSteps, unsigned nCores);

    unsigned dataRows() const;
    unsigned dataCols() const;

    float chi2() const;
    uint64_t nAtoms() const;

    #ifdef GAPS_DEBUG
    float getAvgQueue() const;
    bool internallyConsistent();
    #endif

    // serialization
    friend Archive& operator<< <T, MatA, MatB> (Archive &ar, GibbsSampler &sampler);
    friend Archive& operator>> <T, MatA, MatB> (Archive &ar, GibbsSampler &sampler);
};

class AmplitudeGibbsSampler : public GibbsSampler<AmplitudeGibbsSampler, ColMatrix, RowMatrix>
{
private:

    friend class GibbsSampler;
    friend class PatternGibbsSampler;

    void updateAPMatrix(unsigned row, unsigned col, float delta);

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2,
        unsigned c2);

    float computeDeltaLL(unsigned row, unsigned col, float mass);
    float computeDeltaLL(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);

public:

    template <class DataType>
    AmplitudeGibbsSampler(const DataType &data, bool transposeData,
        unsigned nPatterns, bool partitionRows,
        const std::vector<unsigned> &indices);

    void sync(PatternGibbsSampler &sampler);
};

class PatternGibbsSampler : public GibbsSampler<PatternGibbsSampler, RowMatrix, ColMatrix>
{
private:

    friend class GibbsSampler;
    friend class AmplitudeGibbsSampler;

    void updateAPMatrix(unsigned row, unsigned col, float delta);

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2,
        unsigned c2);

    float computeDeltaLL(unsigned row, unsigned col, float mass);
    float computeDeltaLL(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);

public:

    template <class DataType>
    PatternGibbsSampler(const DataType &data, bool transposeData,
        unsigned nPatterns, bool partitionRows,
        const std::vector<unsigned> &indices);

    void sync(AmplitudeGibbsSampler &sampler);
};

/******************** IMPLEMENTATION OF TEMPLATED FUNCTIONS *******************/

template <class DataType>
AmplitudeGibbsSampler::AmplitudeGibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool partitionRows,
const std::vector<unsigned> &indices)
    :
GibbsSampler(data, transposeData, nPatterns, true, partitionRows, indices)
{}

template <class DataType>
PatternGibbsSampler::PatternGibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool partitionRows,
const std::vector<unsigned> &indices)
    :
GibbsSampler(data, transposeData, nPatterns, false, partitionRows, indices)
{}

template <class T, class MatA, class MatB>
template <class DataType>
GibbsSampler<T, MatA, MatB>::GibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool amp, bool partitionRows,
const std::vector<unsigned> &indices)
    :
mDMatrix(data, transposeData, partitionRows, indices),
mSMatrix(mDMatrix.pmax(0.1f, 0.1f)),
mAPMatrix(mDMatrix.nRow(), mDMatrix.nCol()),
mMatrix(amp ? mDMatrix.nRow() : mDMatrix.nCol(), nPatterns),
mOtherMatrix(NULL),
mDomain(mMatrix.nRow() * mMatrix.nCol()),
mLambda(0.f),
mMaxGibbsMass(100.f),
mAnnealingTemp(1.f),
mNumRows(mMatrix.nRow()),
mNumCols(mMatrix.nCol()),
mAvgQueue(0.f),
mNumQueues(0.f),
mNumBins(mMatrix.nRow() * mMatrix.nCol()),
mBinSize(std::numeric_limits<uint64_t>::max() / mNumBins),
mDomainLength(mBinSize * mNumBins)
{
    // default sparsity parameters
    setSparsity(0.01, false);
}

template <class T, class MatA, class MatB>
template <class DataType>
void GibbsSampler<T, MatA, MatB>::setUncertainty(const DataType &unc,
bool transpose, bool partitionRows, const std::vector<unsigned> &indices)
{
    mSMatrix = MatB(unc, transpose, partitionRows, indices);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setSparsity(float alpha, bool singleCell)
{
    float meanD = singleCell ? gaps::algo::nonZeroMean(mDMatrix) :
        gaps::algo::mean(mDMatrix);

    unsigned nPatterns = mMatrix.nCol();

    mAlpha = alpha;
    mLambda = alpha * std::sqrt(nPatterns / meanD);
}

template <class T, class MatA, class MatB>
unsigned GibbsSampler<T, MatA, MatB>::getRow(uint64_t pos) const
{
    return pos / (mBinSize * mNumCols);
}

template <class T, class MatA, class MatB>
unsigned GibbsSampler<T, MatA, MatB>::getCol(uint64_t pos) const
{
    return (pos / mBinSize) % mNumCols;
}

template <class T, class MatA, class MatB>
bool GibbsSampler<T, MatA, MatB>::canUseGibbs(unsigned col) const
{
    return !gaps::algo::isVectorZero(mOtherMatrix->colPtr(col),
        mOtherMatrix->nRow());
}

template <class T, class MatA, class MatB>
bool GibbsSampler<T, MatA, MatB>::canUseGibbs(unsigned c1, unsigned c2) const
{
    return canUseGibbs(c1) || canUseGibbs(c2);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setMaxGibbsMass(float max)
{
    mMaxGibbsMass = max;
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setMatrix(const Matrix &mat)
{   
    GAPS_ASSERT(mat.nRow() == mMatrix.nRow());
    GAPS_ASSERT(mat.nCol() == mMatrix.nCol());
    mMatrix = mat;
}

template <class T, class MatA, class MatB>
unsigned GibbsSampler<T, MatA, MatB>::dataRows() const
{
    return mDMatrix.nRow();
}

template <class T, class MatA, class MatB>
unsigned GibbsSampler<T, MatA, MatB>::dataCols() const
{
    return mDMatrix.nCol();
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

template <class T, class MatA, class MatB>
T* GibbsSampler<T, MatA, MatB>::impl()
{
    return static_cast<T*>(this);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::update(unsigned nSteps, unsigned nCores)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        makeAndProcessProposal();
        ++n;
    }
}

template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::deathProb(uint64_t nAtoms) const
{
    double size = static_cast<double>(mDomainLength);
    double term1 = (size - static_cast<double>(nAtoms)) / size;
    double term2 = mAlpha * static_cast<double>(mNumBins) * term1;
    return static_cast<double>(nAtoms) / (static_cast<double>(nAtoms) + term2);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::makeAndProcessProposal()
{
    // always birth when no atoms exist
    if (mDomain.size() == 0)
    {
        return birth();
    }

    float bdProb = mDomain.size() < 2 ? 0.6667f : 0.5f;

    float u1 = mPropRng.uniform();
    float u2 = mPropRng.uniform();

    float lowerBound = deathProb(mDomain.size());

    if (u1 <= bdProb)
    {
        return u2 < lowerBound ? death() : birth();
    }
    return (u1 < 0.75f || mDomain.size() < 2) ? move() : exchange();
}

// add an atom at pos, calculate mass either with an exponential distribution
// or with the gibbs mass calculation
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::birth()
{
    uint64_t pos = mDomain.randomFreePosition();
    AtomicProposal prop = AtomicProposal('B', pos);

    unsigned row = getRow(prop.birthPos);
    unsigned col = getCol(prop.birthPos);

    // calculate proposed mass
    float mass = 0.f;
    if (canUseGibbs(col))
    {
        AlphaParameters alpha = impl()->alphaParameters(row, col);
        mass = gibbsMass(alpha, &(prop.rng)).value; // 0 if it fails
    }
    else
    {
        mass = prop.rng.exponential(mLambda);
    }

    // accept mass as long as it's non-zero
    if (mass >= gaps::epsilon)
    {
        mDomain.insert(prop.birthPos, mass);
        mMatrix(row, col) += mass;
        impl()->updateAPMatrix(row, col, mass);
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass calculation
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::death()
{
    Atom* a = mDomain.randomAtom();
    AtomicProposal prop = AtomicProposal('D', a);

    // calculate bin for this atom
    unsigned row = getRow(prop.atom1->pos);
    unsigned col = getCol(prop.atom1->pos);

    // kill off atom
    float newVal = gaps::max(mMatrix(row, col) - prop.atom1->mass, 0.f);
    impl()->updateAPMatrix(row, col, newVal - mMatrix(row, col));
    mMatrix(row, col) = newVal;

    // calculate rebirth mass
    float rebirthMass = prop.atom1->mass;
    AlphaParameters alpha = impl()->alphaParameters(row, col);
    if (canUseGibbs(col))
    {
        OptionalFloat gMass = gibbsMass(alpha, &(prop.rng));
        if (gMass.hasValue)
        {
            rebirthMass = gMass.value;
        }
    }

    // accept/reject rebirth
    float deltaLL = rebirthMass * (alpha.su - alpha.s * rebirthMass / 2.f);
    if (deltaLL * mAnnealingTemp >= std::log(prop.rng.uniform()))
    {
        prop.atom1->mass = rebirthMass;
        mMatrix(row, col) += rebirthMass;
        impl()->updateAPMatrix(row, col, rebirthMass);
    }
    else
    {
        mDomain.erase(prop.atom1->pos);
    }
}

// move mass from src to dest in the atomic domain
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::move()
{
    AtomNeighborhood hood = mDomain.randomAtomWithNeighbors();
    uint64_t lbound = hood.hasLeft() ? hood.left->pos : 0;
    uint64_t rbound = hood.hasRight() ? hood.right->pos : mDomainLength;
    uint64_t newLocation = mPropRng.uniform64(lbound + 1, rbound - 1);

    unsigned r1 = getRow(hood.center->pos);
    unsigned c1 = getCol(hood.center->pos);
    unsigned r2 = getRow(newLocation);
    unsigned c2 = getCol(newLocation);

    if (r1 != r2 || c1 != c2)
    {
        AtomicProposal prop = AtomicProposal('M', hood.center, newLocation);

        float deltaLL = impl()->computeDeltaLL(r1, c1, -1.f * prop.atom1->mass,
            r2, c2, prop.atom1->mass);
        if (deltaLL * mAnnealingTemp > std::log(prop.rng.uniform()))
        {
            prop.atom1->pos = prop.moveDest;

            float newVal = gaps::max(mMatrix(r1, c1) - prop.atom1->mass, 0.f);
            impl()->updateAPMatrix(r1, c1, newVal - mMatrix(r1, c1));
            mMatrix(r1, c1) = newVal;

            mMatrix(r2, c2) += prop.atom1->mass;
            impl()->updateAPMatrix(r2, c2, prop.atom1->mass);
        }
    }
    else
    {
        hood.center->pos = newLocation;
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small after
// the exchange
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::exchange()
{
    AtomNeighborhood hood = mDomain.randomAtomWithRightNeighbor();
    Atom* a1 = hood.center;
    Atom* a2 = hood.hasRight() ? hood.right : mDomain.front();

    unsigned r1 = getRow(a1->pos);
    unsigned c1 = getCol(a1->pos);
    unsigned r2 = getRow(a2->pos);
    unsigned c2 = getCol(a2->pos);

    if (r1 != r2 || c1 != c2)
    {
        AtomicProposal prop = AtomicProposal('E', a1, a2);

        float m1 = prop.atom1->mass;
        float m2 = prop.atom2->mass;

        if (canUseGibbs(c1, c2))
        {
            AlphaParameters alpha = impl()->alphaParameters(r1, c1, r2, c2);
            OptionalFloat gMass = gibbsMass(alpha, m1, m2, &(prop.rng));
            if (gMass.hasValue)
            {
                acceptExchange(&prop, gMass.value, r1, c1, r2, c2);
                return;
            }
        }

        float newMass = prop.rng.truncGammaUpper(m1 + m2, 2.f, 1.f / mLambda);

        float delta = m1 > m2 ? newMass - m1 : m2 - newMass; // change larger mass
        float pOldMass = 2.f * newMass > m1 + m2 ? gaps::max(m1, m2) : gaps::min(m1, m2);

        float pNew = gaps::d_gamma(newMass, 2.f, 1.f / mLambda);
        float pOld = gaps::d_gamma(pOldMass, 2.f, 1.f / mLambda);

        if (pOld == 0.f && pNew != 0.f) // special case
        {
            acceptExchange(&prop, delta, r1, c1, r2, c2);
            return;
        }

        float deltaLL = impl()->computeDeltaLL(r1, c1, delta, r2, c2, -delta);
        float priorLL = (pOld == 0.f) ? 1.f : pOld / pNew;
        float u = std::log(prop.rng.uniform() * priorLL);
        if (u < deltaLL * mAnnealingTemp)
        {
            acceptExchange(&prop, delta, r1, c1, r2, c2);
            return;
        }
    }
}

// helper function for acceptExchange, used to conditionally update the mass
// at a single position, deleting the atom if the new mass is too small
template <class T, class MatA, class MatB>
bool GibbsSampler<T, MatA, MatB>::updateAtomMass(Atom *atom, float delta)
{
    if (atom->mass + delta < gaps::epsilon)
    {
        DEBUG_PING // want to know if this ever happens
        mDomain.erase(atom->pos);
        return false;
    }
    atom->mass += delta;
    return true;
}

// helper function for exchange step, updates the atomic domain, matrix, and
// A*P matrix
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::acceptExchange(AtomicProposal *prop,
float d1, unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    float d2 = -1.f * d1;
    bool b1 = updateAtomMass(prop->atom1, d1);
    bool b2 = updateAtomMass(prop->atom2, d2);
    GAPS_ASSERT(b1 || b2);
    
    // delete entire atom if resize would make it too small
    if (!b1) { d1 = -1.f * prop->atom1->mass; }
    if (!b2) { d2 = -1.f * prop->atom2->mass; }

    // ensure matrix values don't go negative (truncation error at fault)
    float v1 = gaps::max(mMatrix(r1, c1) + d1, 0.f);
    impl()->updateAPMatrix(r1, c1, v1 - mMatrix(r1, c1));
    mMatrix(r1, c1) = v1;


    float v2 = gaps::max(mMatrix(r2, c2) + d2, 0.f);
    impl()->updateAPMatrix(r2, c2, v2 - mMatrix(r2, c2));
    mMatrix(r2, c2) = v2;
}

template <class T, class MatA, class MatB>
OptionalFloat GibbsSampler<T, MatA, MatB>::gibbsMass(AlphaParameters alpha,
GapsRng *rng)
{        
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = (alpha.su - mLambda) / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::p_norm(0.f, mean, sd);

        if (pLower < 1.f)
        {
            float m = rng->inverseNormSample(pLower, 1.f, mean, sd);
            float gMass = gaps::min(m, mMaxGibbsMass / mLambda);
            if (gMass >= gaps::epsilon)
            {
                return OptionalFloat(gMass);
            }
        }
    }
    return OptionalFloat();
}

template <class T, class MatA, class MatB>
OptionalFloat GibbsSampler<T, MatA, MatB>::gibbsMass(AlphaParameters alpha,
float m1, float m2, GapsRng *rng)
{
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = alpha.su / alpha.s; // lambda cancels out
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::p_norm(-m1, mean, sd);
        float pUpper = gaps::p_norm(m2, mean, sd);

        if (!(pLower >  0.95f || pUpper < 0.05f))
        {
            float delta = rng->inverseNormSample(pLower, pUpper, mean, sd);
            float gMass = gaps::min(gaps::max(-m1, delta), m2); // conserve mass
            return OptionalFloat(gMass);
        }
    }
    return OptionalFloat();
}

#ifdef GAPS_DEBUG
template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::getAvgQueue() const
{
    return mAvgQueue;
}

template <class T, class MatA, class MatB>
bool GibbsSampler<T, MatA, MatB>::internallyConsistent()
{
    return true;
}
#endif // GAPS_DEBUG

template <class T, class MatA, class MatB>
Archive& operator<<(Archive &ar, GibbsSampler<T, MatA, MatB> &s)
{
    ar << s.mMatrix << s.mAPMatrix << s.mDomain << s.mLambda
        << s.mMaxGibbsMass << s.mAnnealingTemp << s.mNumRows << s.mNumCols
        << s.mBinSize << s.mAvgQueue << s.mNumQueues;
    return ar;
}

template <class T, class MatA, class MatB>
Archive& operator>>(Archive &ar, GibbsSampler<T, MatA, MatB> &s)
{
    ar >> s.mMatrix >> s.mAPMatrix >> s.mDomain >> s.mLambda
        >> s.mMaxGibbsMass >> s.mAnnealingTemp >> s.mNumRows >> s.mNumCols
        >> s.mBinSize >> s.mAvgQueue >> s.mNumQueues;
    return ar;
}

#endif // __COGAPS_GIBBS_SAMPLER_H__