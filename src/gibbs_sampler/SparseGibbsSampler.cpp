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
    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));

    TemplatedSparseIterator<2> it(mDMatrix.getCol(row), mOtherMatrix->getCol(col));
    while (!it.atEnd())
    {
        s += (get<2>(it) * get<2>(it)) / (get<1>(it) * get<1>(it))
            - (get<2>(it) * get<2>(it));

        s_mu += get<2>(it) / get<1>(it)
            + (get<2>(it) - get<2>(it) / (get<1>(it) * get<1>(it)))
            * gaps::dot(mMatrix.getRow(row), mOtherMatrix->getRow(it.getIndex()));
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
        float s = -2.f * mBeta * mZ2(c1,c2) + a1.s + a2.s;

        TemplatedSparseIterator<3> it(mDMatrix.getCol(r1), mOtherMatrix->getCol(c1),
            mOtherMatrix->getCol(c2));
        while (!it.atEnd())
        {
            float term1 = 2.f * get<2>(it) * get<3>(it);
            s += mBeta * (term1 - term1 / GAPS_SQ(get<1>(it)));
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
    float s = mZ1[col];
    float s_mu = -1.f * gaps::dot(mMatrix.getRow(row), mZ2.getCol(col));
    s_mu -= ch * mZ2(col,col);

    TemplatedSparseIterator<2> it(mDMatrix.getCol(row), mOtherMatrix->getCol(col));
    while (!it.atEnd())
    {
        s += (get<2>(it) * get<2>(it))
            / (get<1>(it) * get<1>(it))
            - (get<2>(it) * get<2>(it));

        s_mu += get<2>(it) / get<1>(it)
            + (get<2>(it) - get<2>(it) / (get<1>(it) * get<1>(it)))
            * (gaps::dot(mMatrix.getRow(row), mOtherMatrix->getRow(it.getIndex()))
            + ch * mOtherMatrix->operator()(it.getIndex(), col));
        it.next();
    }
    return AlphaParameters(s, s_mu) * mBeta;
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
