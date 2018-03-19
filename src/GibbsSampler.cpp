#include "GibbsSampler.h"

AmplitudeGibbsSampler::AmplitudeGibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha, float maxGibbsMass)
    : GibbsSampler(D, S, D.nrow(), nFactor, alpha)
{
    float meanD = gaps::algo::mean(mDMatrix);
    mLambda = alpha * std::sqrt(nFactor / meanD);
    mMaxGibbsMass = maxGibbsMass / mLambda;
}

unsigned AmplitudeGibbsSampler::getRow(uint64_t pos) const
{
    unsigned row = pos / (mBinSize * mNumCols);
    GAPS_ASSERT(row >= 0 && row < mNumRows);
    return row;
}

unsigned AmplitudeGibbsSampler::getCol(uint64_t pos) const
{
    unsigned col = (pos / mBinSize) % mNumCols;
    GAPS_ASSERT(col >= 0 && col < mNumCols);
    return col;
}

bool AmplitudeGibbsSampler::canUseGibbs(unsigned row, unsigned col) const
{
    return !gaps::algo::isRowZero(*mOtherMatrix, col);
}

bool AmplitudeGibbsSampler::canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const
{
    return !gaps::algo::isRowZero(*mOtherMatrix, c1)
        && !gaps::algo::isRowZero(*mOtherMatrix, c2);
}

void AmplitudeGibbsSampler::sync(PatternGibbsSampler &sampler)
{
    mOtherMatrix = &sampler.mMatrix;
    mAPMatrix = sampler.mAPMatrix;
}

void AmplitudeGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    for (unsigned j = 0; j < mAPMatrix.nCol(); ++j)
    {
        mAPMatrix(row,j) += delta * (*mOtherMatrix)(j,col);
    }
}

PatternGibbsSampler::PatternGibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha, float maxGibbsMass)
    : GibbsSampler(D, S, nFactor, D.ncol(), alpha)
{
    float meanD = gaps::algo::mean(mDMatrix);
    mLambda = alpha * std::sqrt(nFactor / meanD);
    mMaxGibbsMass = maxGibbsMass / mLambda;
}

unsigned PatternGibbsSampler::getRow(uint64_t pos) const
{
    unsigned row = (pos / mBinSize) % mNumRows;
    GAPS_ASSERT(row >= 0 && row < mNumRows);
    return row;
}

unsigned PatternGibbsSampler::getCol(uint64_t pos) const
{
    unsigned col = pos / (mBinSize * mNumRows);
    GAPS_ASSERT(col >= 0 && col < mNumCols);
    return col;
}

bool PatternGibbsSampler::canUseGibbs(unsigned row, unsigned col) const
{
    return !gaps::algo::isColZero(*mOtherMatrix, row);
}

bool PatternGibbsSampler::canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const
{
    return !gaps::algo::isColZero(*mOtherMatrix, r1)
        && !gaps::algo::isColZero(*mOtherMatrix, r2);
}

void PatternGibbsSampler::sync(AmplitudeGibbsSampler &sampler)
{
    mOtherMatrix = &sampler.mMatrix;
    mAPMatrix = sampler.mAPMatrix;
}

void PatternGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    for (unsigned i = 0; i < mAPMatrix.nRow(); ++i)
    {
        mAPMatrix(i,col) += delta * (*mOtherMatrix)(i,row);
    }
}