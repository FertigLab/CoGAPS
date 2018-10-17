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
        Vector S(gaps::pmax(D, 0.1f));
        for (unsigned i = 0; i < D.size(); ++i)
        {
            float ap = gaps::dot(mMatrix.getRow(j), mOtherMatrix->getRow(i));
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
    float s = Z1(col);
    float s_mu = 0.f;
    for (unsigned i = 0; i < mNumPatterns; ++i)
    {
        s_mu += mMatrix(row,i) * Z2(col,i);
    }
    s_mu *= -1.f;

    SparseIteratorTwo it(mDMatrix.getCol(row), mOtherMatrix->getCol(col));
    while (!it.atEnd())
    {
        float term1 = it.getValue_2() / it.getValue_1();
        float term2 = it.getValue_2() - term1 / it.getValue_1();
        s += term1 * term1 - it.getValue_2() * it.getValue_2();
        s_mu += term1 + term2 * gaps::dot(mMatrix.getRow(row),
            mOtherMatrix->getRow(it.getIndex()));
        it.next();
    }
    return AlphaParameters(s, s_mu) * mBeta;
}

AlphaParameters SparseGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        AlphaParameters a1 = alphaParameters(r1, c1);
        AlphaParameters a2 = alphaParameters(r2, c2);
        float s = -2.f * mBeta * Z2(c1,c2) + a1.s + a2.s;

        SparseIteratorThree it(mDMatrix.getCol(r1), mOtherMatrix->getCol(c1),
            mOtherMatrix->getCol(c2));
        while (!it.atEnd())
        {
            float term1 = 2.f * it.getValue_2() * it.getValue_3();
            s += mBeta * (term1 - term1 / GAPS_SQ(it.getValue_1()));
            it.next();
        }
        GAPS_ASSERT(s >= 0.f);
        return AlphaParameters(s, a1.s_mu - a2.s_mu);
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

AlphaParameters SparseGibbsSampler::alphaParametersWithChange(unsigned row,
unsigned col, float ch)
{
    float s = Z1(col);
    float s_mu = 0.f;
    for (unsigned i = 0; i < mNumPatterns; ++i)
    {
       s_mu += mMatrix(row,i) * Z2(col,i);
    }
    s_mu += ch * Z2(col,col);
    s_mu *= -1.f;

    SparseIteratorTwo it(mDMatrix.getCol(row), mOtherMatrix->getCol(col));
    while (!it.atEnd())
    {
       float term1 = it.getValue_2() / it.getValue_1();
       float term2 = it.getValue_2() - term1 / it.getValue_1();
       s += term1 * term1 - it.getValue_2() * it.getValue_2();
       s_mu += term1 + term2 * gaps::dot(mMatrix.getRow(row),
           mOtherMatrix->getRow(it.getIndex()));
       s_mu += term2 * mOtherMatrix->operator()(it.getIndex(), col) * ch;
       it.next();
    }
    return AlphaParameters(s, s_mu) * mBeta;
}

float& SparseGibbsSampler::Z1(unsigned pattern)
{
    return mZ1[pattern];
}

float& SparseGibbsSampler::Z2(unsigned pattern1, unsigned pattern2)
{
    if (pattern1 > pattern2)
    {
        unsigned temp = pattern2;
        pattern2 = pattern1;
        pattern1 = temp;
    }
    unsigned offset = pattern1 * (2 * mNumPatterns - pattern1 - 1);
    offset /= 2;
    return mZ2[offset + pattern2];
}

void SparseGibbsSampler::generateLookupTables()
{
    for (unsigned i = 0; i < mNumPatterns; ++i)
    {
        Z1(i) = 0.f;
        for (unsigned k = 0; k < mOtherMatrix->nRow(); ++k)
        {
            Z1(i) += GAPS_SQ(mOtherMatrix->operator()(k,i));
        }
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