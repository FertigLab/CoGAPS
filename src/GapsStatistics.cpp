#include "GapsStatistics.h"

GapsStatistics::GapsStatistics(unsigned nRow, unsigned nCol, unsigned nFactor)
    : mAMeanMatrix(nRow, nFactor), mAStdMatrix(nRow, nFactor),
        mPMeanMatrix(nFactor, nCol), mPStdMatrix(nFactor, nCol), mStatUpdates(0)
{}

void GapsStatistics::update(const AmplitudeGibbsSampler &ASampler,
const PatternGibbsSampler &PSampler)
{
    mStatUpdates++;

    mAMeanMatrix += ASampler.mAMatrix;
    mPMeanMatrix += PSampler.mPMatrix;
    
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
