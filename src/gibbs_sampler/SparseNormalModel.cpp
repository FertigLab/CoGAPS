#include "SparseNormalModel.h"
#include "../data_structures/SparseIterator.h"
#include "../math/Math.h"
#include "../math/Random.h"
#include "../math/MatrixMath.h"
#include "../math/VectorMath.h"
#include "../utils/Archive.h"
#include "../utils/GapsAssert.h"

#define GAPS_SQ(x) ((x) * (x))

#define COUNT_LOWER_BITS(u, pos) __builtin_popcountll((u) & ((1ull << (pos)) - 1ull))
#define CLEAR_LOWER_BITS(u, pos) (((pos) == 63) ? 0 : (u) & ~((1ull << ((pos) + 1ull)) - 1ull))
#define COUNT_BITS(u) __builtin_popcountll(u)
#define GET_FIRST_SET_BIT(u) (__builtin_ffsll(u) - 1)

void SparseNormalModel::setMatrix(const Matrix &mat)
{
    mMatrix = mat;
}

void SparseNormalModel::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

void SparseNormalModel::sync(const SparseNormalModel &model, unsigned nThreads) // NOLINT
{
    mOtherMatrix = &(model.mMatrix);
    generateLookupTables();
}

// required for GibbsSampler interface
void SparseNormalModel::extraInitialization()
{
    // nop - not needed
}

float SparseNormalModel::chiSq() const
{
    float chisq = 0.f;
    for (unsigned j = 0; j < mDMatrix.nCol(); ++j)
    {
        for (unsigned i = 0; i < mDMatrix.nRow(); ++i)
        {
            float dot = gaps::dot(mMatrix.getRow(j), mOtherMatrix->getRow(i));
            chisq += dot * dot;
        }

        SparseIterator<1> it(mDMatrix.getCol(j));
        while (!it.atEnd())
        {
            float dot = gaps::dot(mMatrix.getRow(j), mOtherMatrix->getRow(it.getIndex()));
            float dsq = get<1>(it) * get<1>(it);
            chisq += 1 + dot * (dot - 2 * get<1>(it) - dsq * dot) / dsq;
            it.next();
        }
    }
    return chisq * mBeta;
}

float SparseNormalModel::dataSparsity() const
{
    return gaps::sparsity(mDMatrix);
}

uint64_t SparseNormalModel::nElements() const
{
    return mMatrix.nRow() * mMatrix.nCol();
}

uint64_t SparseNormalModel::nPatterns() const
{
    return mMatrix.nCol();
}

float SparseNormalModel::annealingTemp() const
{
    return mAnnealingTemp;
}

float SparseNormalModel::lambda() const
{
    return mLambda;
}

float SparseNormalModel::maxGibbsMass() const
{
    return mMaxGibbsMass;
}

bool SparseNormalModel::canUseGibbs(unsigned col) const
{
    return !gaps::isVectorZero(mOtherMatrix->getCol(col));
}

bool SparseNormalModel::canUseGibbs(unsigned c1, unsigned c2) const
{
    return canUseGibbs(c1) || canUseGibbs(c2);
}

void SparseNormalModel::changeMatrix(unsigned row, unsigned col, float delta)
{
    mMatrix.add(row, col, delta);
    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
}

void SparseNormalModel::safelyChangeMatrix(unsigned row, unsigned col, float delta)
{
    float newVal = gaps::max(mMatrix(row, col) + delta, 0.f);
    mMatrix.set(row, col, newVal);
    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
}

float SparseNormalModel::deltaLogLikelihood(unsigned r1, unsigned c1, unsigned r2,
unsigned c2, float mass)
{
    AlphaParameters alpha = alphaParameters(r1, c1, r2, c2);
    alpha *= mAnnealingTemp;
    return -1.f * mass * (alpha.s_mu + alpha.s * mass / 2.f);
}

OptionalFloat SparseNormalModel::sampleBirth(unsigned row, unsigned col, GapsRng *rng)
{
    AlphaParameters alpha = alphaParameters(row, col);
    return gibbsMass(alpha * mAnnealingTemp, 0.f, mMaxGibbsMass, rng, mLambda);
}

OptionalFloat SparseNormalModel::sampleDeathAndRebirth(unsigned row, unsigned col,
float delta, GapsRng *rng)
{
    AlphaParameters alpha = alphaParametersWithChange(row, col, delta);
    OptionalFloat mass = gibbsMass(alpha * mAnnealingTemp, 0.f, mMaxGibbsMass, rng, mLambda);
    if (mass.hasValue())
    {
        float deltaLL = mass.value() * (alpha.s_mu - alpha.s * mass.value() / 2.f);
        if (std::log(rng->uniform()) < deltaLL)
        {
            return mass;
        }
    }
    return OptionalFloat();
}

OptionalFloat SparseNormalModel::sampleExchange(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2, GapsRng *rng)
{
    AlphaParameters alpha = alphaParameters(r1, c1, r2, c2);
    return gibbsMass(alpha * mAnnealingTemp, -m1, m2, rng);
}

// PERFORMANCE_CRITICAL
AlphaParameters SparseNormalModel::alphaParameters(unsigned row, unsigned col)
{
    const SparseVector &D(mDMatrix.getCol(row));
    const HybridVector &V(mOtherMatrix->getCol(col));
    const std::vector<uint64_t> &bitflags_D(D.getBitFlags());
    const std::vector<uint64_t> &bitflags_V(V.getBitFlags());
    const std::vector<float> &data(D.getData());

    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));

    unsigned sparseIndex = 0;
    unsigned sz = bitflags_D.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        uint64_t d_flags = bitflags_D[i];
        uint64_t common = d_flags & bitflags_V[i];
        while (common != 0u)
        {
            // find common non-zero index, clear out skipped indices of data
            unsigned index = GET_FIRST_SET_BIT(common);
            sparseIndex += COUNT_LOWER_BITS(d_flags, index);
            d_flags = CLEAR_LOWER_BITS(d_flags, index);
            common &= d_flags;

            // get the needed data
            unsigned v_ndx = 64 * i + index;
            float v_val = V[v_ndx];
            float d_val = data[sparseIndex++];

            // compute terms for s and s_mu
            float term1 = v_val / d_val;
            float term2 = v_val - term1 / d_val;
            s += term1 * term1 - v_val * v_val;
            s_mu += term1 + term2 * gaps::dot(mMatrix.getRow(row),
                mOtherMatrix->getRow(v_ndx));
        }
        sparseIndex += COUNT_BITS(d_flags); // skip over any remaining indices
    }
    return AlphaParameters(s, s_mu) * mBeta;
}

// PERFORMANCE_CRITICAL
AlphaParameters SparseNormalModel::alphaParametersWithChange(unsigned row,
unsigned col, float ch)
{
    const SparseVector &D(mDMatrix.getCol(row));
    const HybridVector &V(mOtherMatrix->getCol(col));
    const std::vector<uint64_t> &bitflags_D(D.getBitFlags());
    const std::vector<uint64_t> &bitflags_V(V.getBitFlags());
    const std::vector<float> &data(D.getData());

    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));
    s_mu -= ch * mZ2(col,col);

    unsigned sparseIndex = 0;
    unsigned sz = bitflags_D.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        uint64_t d_flags = bitflags_D[i];
        uint64_t common = d_flags & bitflags_V[i];
        while (common != 0u)
        {
            // find common non-zero index, clear out skipped indices of data
            unsigned index = GET_FIRST_SET_BIT(common);
            sparseIndex += COUNT_LOWER_BITS(d_flags, index);
            d_flags = CLEAR_LOWER_BITS(d_flags, index);
            common &= d_flags;

            // get the needed data
            unsigned v_ndx = 64 * i + index;
            float v_val = V[v_ndx];
            float d_val = data[sparseIndex++];

            // compute terms for s and s_mu
            float term1 = v_val / d_val;
            float term2 = v_val - term1 / d_val;
            s += term1 * term1 - v_val * v_val;
            s_mu += term1 + term2 * gaps::dot(mMatrix.getRow(row),
                mOtherMatrix->getRow(v_ndx));
            s_mu += term2 * mOtherMatrix->operator()(v_ndx, col) * ch;
        }
        sparseIndex += COUNT_BITS(d_flags);
    }
    return AlphaParameters(s, s_mu) * mBeta;
}

// PERFORMANCE_CRITICAL
AlphaParameters SparseNormalModel::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        const SparseVector &D(mDMatrix.getCol(r1));
        const HybridVector &V1(mOtherMatrix->getCol(c1));
        const HybridVector &V2(mOtherMatrix->getCol(c2));
        const std::vector<uint64_t> &bitflags_D(D.getBitFlags());
        const std::vector<uint64_t> &bitflags_V1(V1.getBitFlags());
        const std::vector<uint64_t> &bitflags_V2(V2.getBitFlags());
        const std::vector<float> &data(D.getData());

        float s = mZ1[c1] - 2.f * mZ2(c1,c2) + mZ1[c2];
        float s_mu = -1.f * gaps::dot_diff(mMatrix.getRow(r1), mZ2.getCol(c1),
            mZ2.getCol(c2));

        unsigned sparseIndex = 0;
        unsigned sz = bitflags_D.size();
        for (unsigned i = 0; i < sz; ++i)
        {
            uint64_t d_flags = bitflags_D[i];
            uint64_t common = d_flags & (bitflags_V1[i] | bitflags_V2[i]);
            while (common != 0u)
            {
                // find common non-zero index, clear out skipped indices of data
                unsigned index = GET_FIRST_SET_BIT(common);
                sparseIndex += COUNT_LOWER_BITS(d_flags, index);
                d_flags = CLEAR_LOWER_BITS(d_flags, index);
                common &= d_flags;

                // get the needed data
                unsigned v_ndx = 64 * i + index;
                float v1_val = V1[v_ndx];
                float v2_val = V2[v_ndx];
                float d_val = data[sparseIndex++];

                float d_recip = 1.f / d_val;
                float term1 = 1.f - d_recip * d_recip;
                float v_diff = v1_val - v2_val;
                float ap = gaps::dot(mMatrix.getRow(r1), mOtherMatrix->getRow(v_ndx));

                s -= v_diff * v_diff * term1;
                s_mu += v_diff * (ap * term1 + d_recip);
            }
            sparseIndex += COUNT_BITS(d_flags);
        }
        return AlphaParameters(s, s_mu) * mBeta;
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

void SparseNormalModel::generateLookupTables()
{
    unsigned nPatterns = mZ1.size();
    for (unsigned i = 0; i < nPatterns; ++i)
    {
        mZ1[i] = 0.f;
        for (unsigned k = 0; k < mOtherMatrix->nRow(); ++k)
        {
            mZ1[i] += GAPS_SQ(mOtherMatrix->operator()(k,i));
        }
        for (unsigned j = i; j < nPatterns; ++j)
        {
            float d = gaps::dot(mOtherMatrix->getCol(i), mOtherMatrix->getCol(j));
            mZ2(i,j) = d;
            mZ2(j,i) = d;
        }
    }
}

Archive& operator<<(Archive &ar, const SparseNormalModel &m)
{
    ar << m.mMatrix << m.mBeta;
    return ar;
}

Archive& operator>>(Archive &ar, SparseNormalModel &m)
{
    ar >> m.mMatrix >> m.mBeta;
    return ar;
}

