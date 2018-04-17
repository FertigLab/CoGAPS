#include "GapsStatistics.h"
#include "Algorithms.h"

GapsStatistics::GapsStatistics(unsigned nRow, unsigned nCol, unsigned nFactor)
    : mAMeanMatrix(nRow, nFactor), mAStdMatrix(nRow, nFactor),
        mPMeanMatrix(nFactor, nCol), mPStdMatrix(nFactor, nCol),
        mNumPatterns(nFactor), mStatUpdates(0)
{}

void GapsStatistics::update(const AmplitudeGibbsSampler &ASampler,
const PatternGibbsSampler &PSampler)
{
    mStatUpdates++;

    Vector normVec(mNumPatterns);
    for (unsigned j = 0; j < mNumPatterns; ++j)
    {
        normVec[j] = gaps::algo::sum(PSampler.mMatrix.getRow(j));
        normVec[j] = normVec[j] == 0 ? 1.f : normVec[j];

        Vector quot(PSampler.mMatrix.getRow(j) / normVec[j]);
        mPMeanMatrix.getRow(j) += quot;
        mPStdMatrix.getRow(j) += gaps::algo::elementSq(quot);

        Vector prod(ASampler.mMatrix.getCol(j) * normVec[j]);
        mAMeanMatrix.getCol(j) += prod;
        mAStdMatrix.getCol(j) += gaps::algo::elementSq(prod); 
    }
}

Rcpp::NumericMatrix GapsStatistics::AMean() const
{
    return (mAMeanMatrix / mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GapsStatistics::AStd() const
{
    return gaps::algo::computeStdDev(mAStdMatrix, mAMeanMatrix,
        mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GapsStatistics::PMean() const
{
    return (mPMeanMatrix / mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GapsStatistics::PStd() const
{
    return gaps::algo::computeStdDev(mPStdMatrix, mPMeanMatrix,
        mStatUpdates).rMatrix();
}

float GapsStatistics::meanChiSq() const
{
    //ColMatrix A = mAMeanMatrix / mStatUpdates;
    //RowMatrix P = mPMeanMatrix / mStatUpdates;
    //RowMatrix M(A * P);
    return 0.f;
}

