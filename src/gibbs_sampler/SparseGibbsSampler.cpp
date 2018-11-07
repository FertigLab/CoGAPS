#include "SparseGibbsSampler.h"
#include "../data_structures/SparseIterator.h"
#include "../math/VectorMath.h"

#include <bitset>

#define GAPS_SQ(x) ((x) * (x))

/////////// SparseGibbsSampler Function Definitions //////////////

float SparseGibbsSampler::chiSq() const
{
    float chisq = 0.f;
    for (unsigned j = 0; j < mDMatrix.nCol(); ++j)
    {
        Vector D(mDMatrix.getCol(j).getDense());
        Vector S(gaps::pmax(D, 0.1f));
        for (unsigned i = 0; i < D.size(); ++i)
        {
            //float ap = gaps::dot(mMatrix.getRow(j), mOtherMatrix->getRow(i));
            float ap = 0.f;
            chisq += GAPS_SQ(D[i] - ap) / GAPS_SQ(S[i]);
        }
    }
    return chisq;
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

//#define __USE_ITERATOR__
//#define __LOOP_OVER_INDICES__
#define __INLINED_ITERATOR__

AlphaParameters SparseGibbsSampler::alphaParameters(unsigned row, unsigned col)
{
#if defined(__USE_ITERATOR__)

    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));

    SparseIterator<2> it(mDMatrix.getCol(row), mOtherMatrix->getCol(col));
    while (!it.atEnd())
    {
        float term1 = get<2>(it) / get<1>(it);
        float term2 = get<2>(it) - term1 / get<1>(it);
        s += term1 * term1 - get<2>(it) * get<2>(it);
        s_mu += term1 + term2 * gaps::dot(mMatrix.getRow(row),
            mOtherMatrix->getRow(it.getIndex()));
        it.next();
    }
    return AlphaParameters(s, s_mu) * mBeta;

#elif defined(__LOOP_OVER_INDICES__)

    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));

    const std::vector<unsigned> &indices(mDMatrix.getCol(row).getIndices());
    const std::vector<float> &data(mDMatrix.getCol(row).getData());
    const HybridVector &vec(mOtherMatrix->getCol(col));
    unsigned sz = indices.size();

    for (unsigned i = 0; i < sz; ++i)
    {
        if (vec[indices[i]] > 0.f)
        {
            float term1 = vec[indices[i]] / data[i];
            float term2 = vec[indices[i]] - term1 / data[i];
            s += term1 * term1 - vec[indices[i]] * vec[indices[i]];
            s_mu += term1 + term2 * gaps::dot(mMatrix.getRow(row),
                mOtherMatrix->getRow(indices[i]));
        }
    }
    return AlphaParameters(s, s_mu) * mBeta;

#elif defined(__INLINED_ITERATOR__)

    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));

    const SparseVector &D(mDMatrix.getCol(row));
    const HybridVector &V(mOtherMatrix->getCol(col));
    const std::vector<uint64_t> &bitflags_D(D.getBitFlags());
    const std::vector<uint64_t> &bitflags_V(V.getBitFlags());
    const std::vector<float> &data(D.getData());

    unsigned sparseIndex = 0;
    unsigned sz = bitflags_D.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        uint64_t d_flags = bitflags_D[i];
        uint64_t common = d_flags & bitflags_V[i];
        while (common)
        {
            unsigned index = __builtin_ffsll(common) - 1;
            sparseIndex += __builtin_popcountll(d_flags & ((1ull << index) - 1ull));
            d_flags = (index == 63) ? 0 : d_flags & ~((1ull << (index + 1ull)) - 1ull);
            common = d_flags & common;

            float v_val = V[64 * i + index];
            float d_val = data[sparseIndex++];

            float term1 = v_val / d_val;
            float term2 = v_val - term1 / d_val;
            s += term1 * term1 - v_val * v_val;
            s_mu += term1 + term2 * gaps::dot(mMatrix.getRow(row),
                mOtherMatrix->getRow(64 * i + index));
        }
        sparseIndex += __builtin_popcountll(d_flags);
    }
    return AlphaParameters(s, s_mu) * mBeta;

#else

    #error "malformed macros"

#endif
}

AlphaParameters SparseGibbsSampler::alphaParametersWithChange(unsigned row,
unsigned col, float ch)
{
#if defined(__USE_ITERATOR__)

    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));
    s_mu -= ch * mZ2(col,col);

    SparseIterator<2> it(mDMatrix.getCol(row), mOtherMatrix->getCol(col));
    while (!it.atEnd())
    {

        float term1 = get<2>(it) / get<1>(it);
        float term2 = get<2>(it) - term1 / get<1>(it);
        s += term1 * term1 - get<2>(it) * get<2>(it);
        s_mu += term1 + term2 * (gaps::dot(mMatrix.getRow(row),
            mOtherMatrix->getRow(it.getIndex()))
            + ch * mOtherMatrix->operator()(it.getIndex(), col));
        it.next();
    }
    return AlphaParameters(s, s_mu) * mBeta;

#elif defined(__LOOP_OVER_INDICES__)

    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));
    s_mu -= ch * mZ2(col,col);

    const std::vector<unsigned> &indices(mDMatrix.getCol(row).getIndices());
    const std::vector<float> &data(mDMatrix.getCol(row).getData());
    const HybridVector &vec(mOtherMatrix->getCol(col));
    unsigned sz = indices.size();

    for (unsigned i = 0; i < sz; ++i)
    {
        if (vec[indices[i]] > 0.f)
        {
            float term1 = vec[indices[i]] / data[i];
            float term2 = vec[indices[i]] - term1 / data[i];
            s += term1 * term1 - vec[indices[i]] * vec[indices[i]];
            s_mu += term1 + term2 * gaps::dot(mMatrix.getRow(row),
                mOtherMatrix->getRow(indices[i]));
            s_mu += term2 * mOtherMatrix->operator()(indices[i], col) * ch;
        }
    }
    return AlphaParameters(s, s_mu) * mBeta;

#elif defined(__INLINED_ITERATOR__)

    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));
    s_mu -= ch * mZ2(col,col);

    const SparseVector &D(mDMatrix.getCol(row));
    const HybridVector &V(mOtherMatrix->getCol(col));
    const std::vector<uint64_t> &bitflags_D(D.getBitFlags());
    const std::vector<uint64_t> &bitflags_V(V.getBitFlags());
    const std::vector<float> &data(D.getData());

    unsigned sparseIndex = 0;
    unsigned sz = bitflags_D.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        uint64_t d_flags = bitflags_D[i];
        uint64_t common = d_flags & bitflags_V[i];
        while (common)
        {
            unsigned index = __builtin_ffsll(common) - 1;
            sparseIndex += __builtin_popcountll(d_flags & ((1ull << index) - 1ull));
            d_flags = (index == 63) ? 0 : d_flags & ~((1ull << (index + 1ull)) - 1ull);
            common = d_flags & common;

            float v_val = V[64 * i + index];
            float d_val = data[sparseIndex++];

            float term1 = v_val / d_val;
            float term2 = v_val - term1 / d_val;
            s += term1 * term1 - v_val * v_val;
            s_mu += term1 + term2 * gaps::dot(mMatrix.getRow(row),
                mOtherMatrix->getRow(64 * i + index));
            s_mu += term2 * mOtherMatrix->operator()(64 * i + index, col) * ch;
        }
        sparseIndex += __builtin_popcountll(d_flags);
    }
    return AlphaParameters(s, s_mu) * mBeta;

#else

    #error "malformed macros"

#endif
}

AlphaParameters SparseGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        float s = -1.f * (mZ1[c1] - 2 * mZ2.operator()(c1,c2) + mZ1[c2]);
        float s_mu = -1.f * gaps::dot_shifted(mMatrix.getRow(r1), mZ2.getCol(c1), -mZ1[c1]);

        SparseIterator<3> it(mDMatrix.getCol(r1), mOtherMatrix->getCol(c1), mOtherMatrix->getCol(c2));
        while (!it.atEnd())
        {
            float term1 = get<2>(it) - get<3>(it);

            s += term1 * term1 + (term1 * term1 / get<1>(it) / get<1>(it));
            s_mu += term1 * (get<1>(it) - gaps::dot(mMatrix.getRow(r1), mOtherMatrix->getRow(it.getIndex()))) / get<1>(it) / get<1>(it);
            s_mu += term1 * gaps::dot(mMatrix.getRow(r1), mOtherMatrix->getRow(it.getIndex()));

            it.next();
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
