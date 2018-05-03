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
    return !(gaps::algo::isRowZero(*mOtherMatrix, c1)
        && gaps::algo::isRowZero(*mOtherMatrix, c2));
}

void AmplitudeGibbsSampler::sync(PatternGibbsSampler &sampler)
{
    mOtherMatrix = &(sampler.mMatrix);
    mAPMatrix = sampler.mAPMatrix;
}

void AmplitudeGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    for (unsigned j = 0; j < mAPMatrix.nCol(); ++j)
    {
        mAPMatrix(row,j) += delta * (*mOtherMatrix)(col,j);
    }
}

AlphaParameters AmplitudeGibbsSampler::alphaParameters(unsigned row, unsigned col)
{
    return gaps::algo::alphaParameters(mDMatrix.nCol(), mDMatrix.rowPtr(row),
        mSMatrix.rowPtr(row), mAPMatrix.rowPtr(row), mOtherMatrix->rowPtr(col));
}

AlphaParameters AmplitudeGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        return gaps::algo::alphaParameters(mDMatrix.nCol(), mDMatrix.rowPtr(r1),
            mSMatrix.rowPtr(r1), mAPMatrix.rowPtr(r1), mOtherMatrix->rowPtr(c1),
            mOtherMatrix->rowPtr(c2));
    }
    else
    {
        AlphaParameters a1 = alphaParameters(r1, c1);
        AlphaParameters a2 = alphaParameters(r2, c2);
        return a1 + a2;
    }
}

float AmplitudeGibbsSampler::computeDeltaLL(unsigned row, unsigned col, float mass)
{
    return gaps::algo::deltaLL(mDMatrix.nCol(), mDMatrix.rowPtr(row),
        mSMatrix.rowPtr(row), mAPMatrix.rowPtr(row), mOtherMatrix->rowPtr(col),
        mass);
}

float AmplitudeGibbsSampler::computeDeltaLL(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2)
{
    if (r1 == r2)
    {
        return gaps::algo::deltaLL(mDMatrix.nCol(), mDMatrix.rowPtr(r1),
            mSMatrix.rowPtr(r1), mAPMatrix.rowPtr(r1), mOtherMatrix->rowPtr(c1),
            m1, mOtherMatrix->rowPtr(c2), m2);
    }
    else
    {
        return computeDeltaLL(r1, c1, m1) + computeDeltaLL(r2, c2, m2);
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
    return !(gaps::algo::isColZero(*mOtherMatrix, r1)
        && gaps::algo::isColZero(*mOtherMatrix, r2));
}

void PatternGibbsSampler::sync(AmplitudeGibbsSampler &sampler)
{
    mOtherMatrix = &(sampler.mMatrix);
    mAPMatrix = sampler.mAPMatrix;
}

void PatternGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    for (unsigned i = 0; i < mAPMatrix.nRow(); ++i)
    {
        mAPMatrix(i,col) += delta * (*mOtherMatrix)(i,row);
    }
}

AlphaParameters PatternGibbsSampler::alphaParameters(unsigned row, unsigned col)
{
    return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(col),
        mSMatrix.colPtr(col), mAPMatrix.colPtr(col), mOtherMatrix->colPtr(row));
}

AlphaParameters PatternGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (c1 == c2)
    {
        return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(c1),
            mSMatrix.colPtr(c1), mAPMatrix.colPtr(c1), mOtherMatrix->colPtr(r1),
            mOtherMatrix->colPtr(r2));
    }
    else
    {
        AlphaParameters a1 = alphaParameters(r1, c1);
        AlphaParameters a2 = alphaParameters(r2, c2);
        return a1 + a2;
    }
}

float PatternGibbsSampler::computeDeltaLL(unsigned row, unsigned col, float mass)
{
    return gaps::algo::deltaLL(mDMatrix.nRow(), mDMatrix.colPtr(col),
        mSMatrix.colPtr(col), mAPMatrix.colPtr(col), mOtherMatrix->colPtr(row),
        mass);
}

float PatternGibbsSampler::computeDeltaLL(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2)
{
    if (c1 == c2)
    {
        return gaps::algo::deltaLL(mDMatrix.nRow(), mDMatrix.colPtr(c1),
            mSMatrix.colPtr(c1), mAPMatrix.colPtr(c1), mOtherMatrix->colPtr(r1),
            m1, mOtherMatrix->colPtr(r2), m2);
    }
    else
    {
        return computeDeltaLL(r1, c1, m1) + computeDeltaLL(r2, c2, m2);
    }
}