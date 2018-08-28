#include "GibbsSampler.h"
#include "math/SIMD.h"

/******************** AmplitudeGibbsSampler Implementation ********************/

void AmplitudeGibbsSampler::sync(PatternGibbsSampler &sampler)
{
    mOtherMatrix = &(sampler.mMatrix);
    mAPMatrix = sampler.mAPMatrix;
}

void PatternGibbsSampler::sync(AmplitudeGibbsSampler &sampler)
{
    mOtherMatrix = &(sampler.mMatrix);
    mAPMatrix = sampler.mAPMatrix;
}

void AmplitudeGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->colPtr(col);
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

void PatternGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->colPtr(col);
    float *ap = mAPMatrix.colPtr(row);
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

AlphaParameters AmplitudeGibbsSampler::alphaParameters(unsigned row, unsigned col)
{
    return gaps::algo::alphaParameters(mDMatrix.nCol(), mDMatrix.rowPtr(row),
        mSMatrix.rowPtr(row), mAPMatrix.rowPtr(row), mOtherMatrix->colPtr(col));
}

AlphaParameters AmplitudeGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        return gaps::algo::alphaParameters(mDMatrix.nCol(), mDMatrix.rowPtr(r1),
            mSMatrix.rowPtr(r1), mAPMatrix.rowPtr(r1), mOtherMatrix->colPtr(c1),
            mOtherMatrix->colPtr(c2));
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

AlphaParameters PatternGibbsSampler::alphaParameters(unsigned row, unsigned col)
{
    return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(row),
        mSMatrix.colPtr(row), mAPMatrix.colPtr(row), mOtherMatrix->colPtr(col));
}

AlphaParameters PatternGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(r1),
            mSMatrix.colPtr(r1), mAPMatrix.colPtr(r1), mOtherMatrix->colPtr(c1),
            mOtherMatrix->colPtr(c2));
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

float AmplitudeGibbsSampler::computeDeltaLL(unsigned row, unsigned col, float mass)
{
    AlphaParameters alpha = impl()->alphaParameters(row, col);
    return mass * (alpha.su - alpha.s * mass / 2.f);
}

float AmplitudeGibbsSampler::computeDeltaLL(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2)
{
    if (r1 == r2)
    {
        return gaps::algo::deltaLL(mDMatrix.nCol(), mDMatrix.rowPtr(r1),
            mSMatrix.rowPtr(r1), mAPMatrix.rowPtr(r1), mOtherMatrix->colPtr(c1),
            m1, mOtherMatrix->colPtr(c2), m2);
    }
    return computeDeltaLL(r1, c1, m1) + computeDeltaLL(r2, c2, m2);
}

float PatternGibbsSampler::computeDeltaLL(unsigned row, unsigned col, float mass)
{
    AlphaParameters alpha = impl()->alphaParameters(row, col);
    return mass * (alpha.su - alpha.s * mass / 2.f);
}

float PatternGibbsSampler::computeDeltaLL(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2)
{
    if (r1 == r2)
    {
        return gaps::algo::deltaLL(mDMatrix.nRow(), mDMatrix.colPtr(r1),
            mSMatrix.colPtr(r1), mAPMatrix.colPtr(r1), mOtherMatrix->colPtr(c1),
            m1, mOtherMatrix->colPtr(c2), m2);
    }
    return computeDeltaLL(r1, c1, m1) + computeDeltaLL(r2, c2, m2);
}