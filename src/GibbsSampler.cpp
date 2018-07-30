#include "GibbsSampler.h"
#include "math/SIMD.h"

/******************** AmplitudeGibbsSampler Implementation ********************/

void AmplitudeGibbsSampler::sync(PatternGibbsSampler &sampler)
{
    mOtherMatrix = &(sampler.mMatrix);
    mAPMatrix = sampler.mAPMatrix;
}

void AmplitudeGibbsSampler::recalculateAPMatrix()
{
    mAPMatrix = gaps::algo::matrixMultiplication(mMatrix, *mOtherMatrix);
}

unsigned AmplitudeGibbsSampler::getRow(uint64_t pos) const
{
    return pos / (mBinSize * mNumCols);
}

unsigned AmplitudeGibbsSampler::getCol(uint64_t pos) const
{
    return (pos / mBinSize) % mNumCols;
}

bool AmplitudeGibbsSampler::canUseGibbs(unsigned row, unsigned col) const // NOLINT(misc-unused-parameters)
{
    return !gaps::algo::isVectorZero(mOtherMatrix->rowPtr(col),
        mOtherMatrix->nCol());
}

bool AmplitudeGibbsSampler::canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const
{
    return canUseGibbs(r1, c1) || canUseGibbs(r2, c2);
}

void AmplitudeGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->rowPtr(col);
    float *ap = mAPMatrix.rowPtr(row);
    unsigned size = mAPMatrix.nCol();

    gaps::simd::packedFloat pOther, pAP;
    gaps::simd::Index i(0);
    gaps::simd::packedFloat pDelta(delta);
    for (; i <= size - i.increment(); ++i)
    {
        pOther.load(other + i);
        pAP.load(ap + i);
        pAP += pDelta * pOther;
        pAP.store(ap + i);
    }

    for (unsigned j = i.value(); j < size; ++j)
    {
        ap[j] += delta * other[j];
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
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
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
    return computeDeltaLL(r1, c1, m1) + computeDeltaLL(r2, c2, m2);
}

/********************* PatternGibbsSampler Implementation *********************/

void PatternGibbsSampler::recalculateAPMatrix()
{
    mAPMatrix = gaps::algo::matrixMultiplication(*mOtherMatrix, mMatrix);
}

void PatternGibbsSampler::sync(AmplitudeGibbsSampler &sampler)
{
    mOtherMatrix = &(sampler.mMatrix);
    mAPMatrix = sampler.mAPMatrix;
}

unsigned PatternGibbsSampler::getRow(uint64_t pos) const
{
    return (pos / mBinSize) % mNumRows;
}

unsigned PatternGibbsSampler::getCol(uint64_t pos) const
{
    return pos / (mBinSize * mNumRows);
}

bool PatternGibbsSampler::canUseGibbs(unsigned row, unsigned col) const // NOLINT(misc-unused-parameters)
{
    return !gaps::algo::isVectorZero(mOtherMatrix->colPtr(row),
        mOtherMatrix->nRow());
}

bool PatternGibbsSampler::canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const
{
    return canUseGibbs(r1, c1) || canUseGibbs(r2, c2);
}

void PatternGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->colPtr(row);
    float *ap = mAPMatrix.colPtr(col);
    unsigned size = mAPMatrix.nRow();

    gaps::simd::packedFloat pOther, pAP;
    gaps::simd::Index i(0);
    gaps::simd::packedFloat pDelta(delta);
    for (; i <= size - i.increment(); ++i)
    {
        pOther.load(other + i);
        pAP.load(ap + i);
        pAP += pDelta * pOther;
        pAP.store(ap + i);
    }

    for (unsigned j = i.value(); j < size; ++j)
    {
        ap[j] += delta * other[j];
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
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
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
    return computeDeltaLL(r1, c1, m1) + computeDeltaLL(r2, c2, m2);
}