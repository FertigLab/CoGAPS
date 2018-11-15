#include "SparseGibbsSampler.h"
#include "../data_structures/SparseIterator.h"
#include "../math/VectorMath.h"

#include <bitset>

#define GAPS_SQ(x) ((x) * (x))

#define COUNT_LOWER_BITS(u, pos) __builtin_popcountll((u) & ((1ull << (pos)) - 1ull))
#define CLEAR_LOWER_BITS(u, pos) (((pos) == 63) ? 0 : (u) & ~((1ull << ((pos) + 1ull)) - 1ull))
#define COUNT_BITS(u) __builtin_popcountll(u)
#define GET_FIRST_SET_BIT(u) (__builtin_ffsll(u) - 1)

/////////// SparseGibbsSampler Function Definitions //////////////

float SparseGibbsSampler::chiSq() const
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

void SparseGibbsSampler::sync(const SparseGibbsSampler &sampler, unsigned nThreads)
{
    mOtherMatrix = &(sampler.mMatrix);
    generateLookupTables();
}

// required for GibbsSampler interface
void SparseGibbsSampler::extraInitialization()
{
    // nop - not needed
}

void SparseGibbsSampler::changeMatrix(unsigned row, unsigned col,
float delta)
{
    mMatrix.add(row, col, delta);

    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
}

void SparseGibbsSampler::safelyChangeMatrix(unsigned row,
unsigned col, float delta)
{
    float newVal = gaps::max(mMatrix(row, col) + delta, 0.f);
    mMatrix.add(row, col, newVal - mMatrix(row, col));

    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
}

AlphaParameters SparseGibbsSampler::alphaParameters(unsigned row, unsigned col)
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
        while (common)
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

AlphaParameters SparseGibbsSampler::alphaParametersWithChange(unsigned row,
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
        while (common)
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

AlphaParameters SparseGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
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
            while (common)
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

void SparseGibbsSampler::generateLookupTables()
{
    for (unsigned i = 0; i < mNumPatterns; ++i)
    {
        mZ1[i] = 0.f;
        for (unsigned k = 0; k < mOtherMatrix->nRow(); ++k)
        {
            mZ1[i] += GAPS_SQ(mOtherMatrix->operator()(k,i));
        }
        for (unsigned j = i; j < mNumPatterns; ++j)
        {
            float d = gaps::dot(mOtherMatrix->getCol(i), mOtherMatrix->getCol(j));
            mZ2(i,j) = d;
            mZ2(j,i) = d;
        }
    }
}

Archive& operator<<(Archive &ar, const SparseGibbsSampler &s)
{
    ar << s.mMatrix << s.mDomain << s.mAlpha << s.mLambda << s.mMaxGibbsMass
        << s.mAnnealingTemp << s.mNumPatterns << s.mNumBins << s.mBinLength
        << s.mBeta;
    return ar;
}

Archive& operator>>(Archive &ar, SparseGibbsSampler &s)
{
    ar >> s.mMatrix >> s.mDomain >> s.mAlpha >> s.mLambda >> s.mMaxGibbsMass
        >> s.mAnnealingTemp >> s.mNumPatterns >> s.mNumBins >> s.mBinLength
        >> s.mBeta;
    return ar;
}

#ifdef GAPS_DEBUG
bool SparseGibbsSampler::internallyConsistent() const
{
    return true;
}
#endif
