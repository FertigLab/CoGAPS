#include "GibbsSampler.h"

/*
AmplitudeGibbsSampler::AmplitudeGibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor)
    :
mMatrix(D.nrow(), nFactor), mDMatrix(D), mSMatrix(S), mAPMatrix(D.nrow(), D.ncol()),
mNumRows(D.nrow()), mNumCols(nFactor)
{}

AmplitudeGibbsSampler::AmplitudeGibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha, float maxGibbsmass)
    :
mMatrix(D.nrow(), nFactor), mDMatrix(D), mSMatrix(S), mAPMatrix(D.nrow(), D.ncol()),
mQueue(alpha, D.ncol(), D.nrow() * nFactor), mMaxGibbsMass(maxGibbsmass),
mAnnealingTemp(0), mNumRows(D.nrow()), mNumCols(nFactor)
{}
*/

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
    return !gaps::algo::isRowZero(*mOtherMatrix, col);
}

bool AmplitudeGibbsSampler::canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const
{
    return !gaps::algo::isRowZero(*mOtherMatrix, c1)
        && !gaps::algo::isRowZero(*mOtherMatrix, c2);
}

void AmplitudeGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    for (unsigned j = 0; j < mAPMatrix.nCol(); ++j)
    {
        mAPMatrix(row,j) += delta * (*mOtherMatrix)(j,col);
    }
}

/*
PatternGibbsSampler::PatternGibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor)
    :
mMatrix(nFactor, D.ncol()), mDMatrix(D), mSMatrix(S), mAPMatrix(D.nrow(), D.ncol()),
mNumRows(nFactor), mNumCols(D.ncol())
{}

PatternGibbsSampler::PatternGibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha, float maxGibbsmass)
    :
mMatrix(nFactor, D.ncol()), mDMatrix(D), mSMatrix(S), mAPMatrix(D.nrow(), D.ncol()),
mQueue(alpha, D.nrow(), D.ncol() * nFactor), mMaxGibbsMass(maxGibbsmass),
mAnnealingTemp(0), mNumRows(nFactor), mNumCols(D.ncol())
{}
*/

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
    return !gaps::algo::isColZero(*mOtherMatrix, row);
}

bool PatternGibbsSampler::canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const
{
    return !gaps::algo::isColZero(*mOtherMatrix, r1)
        && !gaps::algo::isColZero(*mOtherMatrix, r2);
}

void PatternGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    for (unsigned i = 0; i < mAPMatrix.nRow(); ++i)
    {
        mAPMatrix(i,col) += delta * (*mOtherMatrix)(i,row);
    }
}
