#include "SparseGibbsSampler.h"
#include "../data_structures/SparseIterator.h"
#include "../math/VectorMath.h"

#define GAPS_SQ(x) ((x) * (x))

/////////// SparseGibbsSampler Function Definitions //////////////

float SparseGibbsSampler::chiSq() const
{
    float chisq = 0.f;
    for (unsigned j = 0; j < mDMatrix.nCol(); ++j)
    {
        Vector D(mDMatrix.getCol(j).getDense());
        for (unsigned i = 0; i < D.size(); ++i)
        {
            float ap = gaps::dot(mMatrix.getRow(j), mOtherMatrix->getRow(i));
            chisq += GAPS_SQ(D[i] - ap) / GAPS_SQ(D[i]);
        }
    }
    return mBeta * chisq;
}

void SparseGibbsSampler::sync(const SparseGibbsSampler &sampler, unsigned nThreads)
{
    mOtherMatrix = &(sampler.mMatrix);
    generateLookupTables();
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

AlphaParameters SparseGibbsSampler::alphaParameters(unsigned row,
unsigned col)
{
    float s = -1.f * Z1(col) * mBeta;
    float s_mu = 0.f;
    for (unsigned i = 0; i < mNumPatterns; ++i)
    {
        s_mu += mMatrix(row,i) * Z2(col,i);
    }

    SparseIteratorTwo it(mDMatrix.getCol(row), mOtherMatrix->getCol(col));
    while (!it.atEnd())
    {
        float term1 = it.getValue_1() / it.getValue_2();
        float term2 = term1 * term1 + it.getValue_1() * it.getValue_1() * mBeta;
        float term3 = mBeta * (it.getValue_1() - term1 / it.getValue_2());
        s += mBeta * term2;
        s_mu += mBeta * term1 + term3 * gaps::dot(mMatrix.getRow(row),
            mOtherMatrix->getRow(it.getIndex()));
        it.next();
    }
    return AlphaParameters(s, s_mu);
}

AlphaParameters SparseGibbsSampler::alphaParameters(unsigned r1,
unsigned c1, unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        AlphaParameters a1 = alphaParameters(r1, c1);
        AlphaParameters a2 = alphaParameters(r2, c2);

        float s = -2.f * mBeta * Z2(c1,c2) + a1.s + a2.s;
        float s_mu = a1.s_mu - a2.s_mu;

        SparseIteratorThree it(mDMatrix.getCol(r1), mOtherMatrix->getCol(c1),
            mOtherMatrix->getCol(c2));
        while (!it.atEnd())
        {
            float term1 = 2.f * it.getValue_1() * it.getValue_2();
            s_mu += mBeta * term1 * (1.f - 1.f / it.getValue_3());
            it.next();
        }
        return AlphaParameters(s, s_mu);
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

AlphaParameters SparseGibbsSampler::alphaParametersWithChange(
unsigned row, unsigned col, float ch)
{
    float s = -1.f * Z1(col) * mBeta;
    float s_mu = 0.f;
    for (unsigned i = 0; i < mNumPatterns; ++i)
    {
        s_mu += mMatrix(row,i) * Z2(col,i);
    }
    s_mu += Z2(col,col) * ch;

    SparseIteratorTwo it(mDMatrix.getCol(row), mOtherMatrix->getCol(col));
    while (!it.atEnd())
    {
        float term1 = it.getValue_1() / it.getValue_2();
        float term2 = term1 * term1 + it.getValue_1() * it.getValue_1() * mBeta;
        float term3 = mBeta * (it.getValue_1() - term1 / it.getValue_2());
        s += mBeta * term2;
        s_mu += mBeta * term1 + term3 * gaps::dot(mMatrix.getRow(row),
            mOtherMatrix->getRow(it.getIndex()));
        s_mu += term3 * mOtherMatrix->operator()(it.getIndex(), col) * ch;
        it.next();
    }
    return AlphaParameters(s, s_mu);
}

float& SparseGibbsSampler::Z1(unsigned pattern)
{
    return mZ1[pattern];
}

float& SparseGibbsSampler::Z2(unsigned pattern1, unsigned pattern2)
{
    unsigned dist = mNumPatterns - pattern1;
    unsigned offset = mNumPatterns * mNumPatterns + mNumPatterns
        - dist * dist + dist;
    offset /= 2;
    return mZ2[offset + pattern2];
}

void SparseGibbsSampler::generateLookupTables()
{
    for (unsigned i = 0; i < mNumPatterns; ++i)
    {
        Z1(i) = gaps::sum(mOtherMatrix->getCol(i));
        for (unsigned j = i; j < mNumPatterns; ++j)
        {
            Z2(i,j) = gaps::dot(mOtherMatrix->getCol(i),
                mOtherMatrix->getCol(j));
        }
    }
}

Archive& operator<<(Archive &ar, SparseGibbsSampler &s)
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