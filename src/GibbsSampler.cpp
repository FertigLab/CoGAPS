#include "GibbsSampler.h"

AmplitudeGibbsSampler::AmplitudeGibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha, float maxGibbsMass)
    :
mAMatrix(D.nrow(), nFactor), mDMatrix(D), mSMatrix(S),
mAPMatrix(D.nrow(), D.ncol()), mQueue(D.nrow() * nFactor, alpha),
mAnnealingTemp(0.f), mNumRows(D.nrow()), mNumCols(nFactor)
{
    mBinSize = std::numeric_limits<uint64_t>::max() / (mNumRows * mNumCols);
    uint64_t remain = std::numeric_limits<uint64_t>::max() % (mNumRows * mNumCols);
    mQueue.setDomainSize(std::numeric_limits<uint64_t>::max() - remain);

    float meanD = gaps::algo::mean(mDMatrix);
    mLambda = alpha * std::sqrt(nFactor / meanD);
    mMaxGibbsMass = maxGibbsMass / mLambda;
}

unsigned AmplitudeGibbsSampler::getRow(uint64_t pos) const
{
    return pos / (mBinSize * mNumCols);
}

unsigned AmplitudeGibbsSampler::getCol(uint64_t pos) const
{
    return (pos / mBinSize) % mNumCols;
}

bool AmplitudeGibbsSampler::canUseGibbs(unsigned row, unsigned col) const
{
    return !gaps::algo::isRowZero(*mPMatrix, col);
}

bool AmplitudeGibbsSampler::canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const
{
    return !gaps::algo::isRowZero(*mPMatrix, c1)
        && !gaps::algo::isRowZero(*mPMatrix, c2);
}

void AmplitudeGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    for (unsigned j = 0; j < mAPMatrix.nCol(); ++j)
    {
        mAPMatrix(row,j) += delta * (*mPMatrix)(j,col);
    }
}

void AmplitudeGibbsSampler::sync(PatternGibbsSampler &sampler)
{
    mPMatrix = &(sampler.mPMatrix);
    mAPMatrix = sampler.mAPMatrix;
}

float AmplitudeGibbsSampler::nAtoms() const
{
    return mDomain.size();
}

void AmplitudeGibbsSampler::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

float AmplitudeGibbsSampler::chi2() const
{
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);   
}

PatternGibbsSampler::PatternGibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha, float maxGibbsMass)
    :
mPMatrix(nFactor, D.ncol()), mDMatrix(D), mSMatrix(S),
mAPMatrix(D.nrow(), D.ncol()), mQueue(nFactor * D.ncol(), alpha),
mAnnealingTemp(0.f), mNumRows(nFactor), mNumCols(D.nrow())
{
    mBinSize = std::numeric_limits<uint64_t>::max() / (mNumRows * mNumCols);
    uint64_t remain = std::numeric_limits<uint64_t>::max() % (mNumRows * mNumCols);
    mQueue.setDomainSize(std::numeric_limits<uint64_t>::max() - remain);

    float meanD = gaps::algo::mean(mDMatrix);
    mLambda = alpha * std::sqrt(nFactor / meanD);
    mMaxGibbsMass = maxGibbsMass / mLambda;
}

unsigned PatternGibbsSampler::getRow(uint64_t pos) const
{
    return (pos / mBinSize) % mNumRows;
}

unsigned PatternGibbsSampler::getCol(uint64_t pos) const
{
    return pos / (mBinSize * mNumRows);
}

bool PatternGibbsSampler::canUseGibbs(unsigned row, unsigned col) const
{
    return !gaps::algo::isColZero(*mAMatrix, row);
}

bool PatternGibbsSampler::canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const
{
    return !gaps::algo::isColZero(*mAMatrix, r1)
        && !gaps::algo::isColZero(*mAMatrix, r2);
}

void PatternGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    for (unsigned i = 0; i < mAPMatrix.nRow(); ++i)
    {
        mAPMatrix(i,col) += delta * (*mAMatrix)(i,row);
    }
}

void PatternGibbsSampler::sync(AmplitudeGibbsSampler &sampler)
{
    mAMatrix = &(sampler.mAMatrix);
    mAPMatrix = sampler.mAPMatrix;
}

float PatternGibbsSampler::nAtoms() const
{
    return mDomain.size();
}

void PatternGibbsSampler::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

float PatternGibbsSampler::chi2() const
{
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);   
}