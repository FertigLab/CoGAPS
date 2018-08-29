#include "GapsStatistics.h"
#include "math/Algorithms.h"

GapsStatistics::GapsStatistics(unsigned nRow, unsigned nCol, unsigned nPatterns)
    :
mAMeanMatrix(nRow, nPatterns), mAStdMatrix(nRow, nPatterns),
mPMeanMatrix(nCol, nPatterns), mPStdMatrix(nCol, nPatterns),
mStatUpdates(0), mNumPatterns(nPatterns)
{}

void GapsStatistics::update(const GibbsSampler &ASampler,
const GibbsSampler &PSampler)
{
    mStatUpdates++;

    // update     
    for (unsigned j = 0; j < mNumPatterns; ++j)
    {
        float norm = gaps::algo::max(PSampler.mMatrix.getCol(j));
        norm = norm == 0.f ? 1.f : norm;

        Vector quot(PSampler.mMatrix.getCol(j) / norm);
        mPMeanMatrix.getCol(j) += quot;
        mPStdMatrix.getCol(j) += gaps::algo::elementSq(quot);

        Vector prod(ASampler.mMatrix.getCol(j) * norm);
        mAMeanMatrix.getCol(j) += prod;
        mAStdMatrix.getCol(j) += gaps::algo::elementSq(prod); // precision loss? use double?
    }
}

ColMatrix GapsStatistics::Amean() const
{
    return mAMeanMatrix / mStatUpdates;
}

ColMatrix GapsStatistics::Asd() const
{
    return gaps::algo::computeStdDev(mAStdMatrix, mAMeanMatrix,
        mStatUpdates);
}

ColMatrix GapsStatistics::Pmean() const
{
    return mPMeanMatrix / mStatUpdates;
}

ColMatrix GapsStatistics::Psd() const
{
    return gaps::algo::computeStdDev(mPStdMatrix, mPMeanMatrix,
        mStatUpdates);
}

float GapsStatistics::meanChiSq(const GibbsSampler &PSampler) const
{
    ColMatrix A = mAMeanMatrix / mStatUpdates;
    ColMatrix P = mPMeanMatrix / mStatUpdates;
    ColMatrix M(gaps::algo::matrixMultiplication(A, P));
    return 2.f * gaps::algo::loglikelihood(PSampler.mDMatrix, PSampler.mSMatrix,
        M);
}

Archive& operator<<(Archive &ar, GapsStatistics &stat)
{
    ar << stat.mAMeanMatrix << stat.mAStdMatrix << stat.mPMeanMatrix
        << stat.mPStdMatrix << stat.mStatUpdates << stat.mNumPatterns;
    return ar;
}

Archive& operator>>(Archive &ar, GapsStatistics &stat)
{
    ar >> stat.mAMeanMatrix >> stat.mAStdMatrix >> stat.mPMeanMatrix
        >> stat.mPStdMatrix >> stat.mStatUpdates >> stat.mNumPatterns;
    return ar;
}

