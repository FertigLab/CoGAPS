#include "DenseGibbsSampler.h"

#define GAPS_SQ(x) ((x) * (x))

float DenseGibbsSampler::chiSq() const
{
    float chisq = 0.f;
    for (unsigned i = 0; i < mDMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < mDMatrix.nCol(); ++j)
        {
            GAPS_ASSERT(mSMatrix(i,j) > 0.f);
            chisq += GAPS_SQ((mDMatrix(i,j) - mAPMatrix(i,j)) / mSMatrix(i,j));
        }
    }
    return chisq;
}

void DenseGibbsSampler::sync(const DenseGibbsSampler &sampler, unsigned nThreads)
{
    // copy transpose of other AP matrix
    GAPS_ASSERT(sampler.mAPMatrix.nRow() == mAPMatrix.nCol());
    GAPS_ASSERT(sampler.mAPMatrix.nCol() == mAPMatrix.nRow());

    unsigned nc = sampler.mAPMatrix.nCol();
    unsigned nr = sampler.mAPMatrix.nRow();

    #pragma omp parallel for num_threads(nThreads)
    for (unsigned j = 0; j < nc; ++j)
    {
        for (unsigned i = 0; i < nr; ++i)
        {
            mAPMatrix(j,i) = sampler.mAPMatrix(i,j);
        }
    }
    mOtherMatrix = &(sampler.mMatrix); // update pointer
}

void DenseGibbsSampler::recalculateAPMatrix()
{
    GAPS_ASSERT(mOtherMatrix->nRow() == mAPMatrix.nRow());
    GAPS_ASSERT(mOtherMatrix->nCol() == mMatrix.nCol());
    GAPS_ASSERT(mMatrix.nRow() == mAPMatrix.nCol());

    for (unsigned j = 0; j < mAPMatrix.nCol(); ++j)
    {
        for (unsigned i = 0; i < mAPMatrix.nRow(); ++i)
        {
            mAPMatrix(i,j) = 0.f;
            for (unsigned k = 0; k < mMatrix.nCol(); ++k)
            {
                mAPMatrix(i,j) += mOtherMatrix->operator()(i,k) * mMatrix(j,k);
            }
        }
    }
}

void DenseGibbsSampler::changeMatrix(unsigned row, unsigned col, float delta)
{
    mMatrix(row, col) += delta;
    updateAPMatrix(row, col, delta);

    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
}

void DenseGibbsSampler::safelyChangeMatrix(unsigned row, unsigned col, float delta)
{
    float newVal = gaps::max(mMatrix(row, col) + delta, 0.f);
    updateAPMatrix(row, col, newVal - mMatrix(row, col));
    mMatrix(row, col) = newVal;
}

// PERFORMANCE_CRITICAL
AlphaParameters DenseGibbsSampler::alphaParameters(unsigned row, unsigned col)
{
    unsigned size = mDMatrix.nRow();
    const float *D = mDMatrix.getCol(row).ptr();
    const float *S = mSMatrix.getCol(row).ptr();
    const float *AP = mAPMatrix.getCol(row).ptr();
    const float *mat = mOtherMatrix->getCol(col).ptr();

    gaps::simd::PackedFloat pMat, pD, pAP, pS;
    gaps::simd::PackedFloat partialS(0.f), partialS_mu(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        gaps::simd::PackedFloat ratio(pMat / (pS * pS));
        partialS += pMat * ratio;
        partialS_mu += ratio * (pD - pAP);
    }

    float s = partialS.scalar(), s_mu = partialS_mu.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        float ratio = mat[j] / (S[j] * S[j]);
        s += mat[j] * ratio;
        s_mu += ratio * (D[j] - AP[j]);
    }
    return AlphaParameters(s, s_mu);
}

// PERFORMANCE_CRITICAL
AlphaParameters DenseGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        unsigned size = mDMatrix.nRow();
        const float *D = mDMatrix.getCol(r1).ptr();
        const float *S = mSMatrix.getCol(r1).ptr();
        const float *AP = mAPMatrix.getCol(r1).ptr();
        const float *mat1 = mOtherMatrix->getCol(c1).ptr();
        const float *mat2 = mOtherMatrix->getCol(c2).ptr();

        gaps::simd::PackedFloat pMat1, pMat2, pD, pAP, pS;
        gaps::simd::PackedFloat partialS(0.f), partialS_mu(0.f);
        gaps::simd::Index i(0);
        for (; i <= size - i.increment(); ++i)
        {   
            pMat1.load(mat1 + i);
            pMat2.load(mat2 + i);
            pD.load(D + i);
            pAP.load(AP + i);
            pS.load(S + i);
            gaps::simd::PackedFloat ratio((pMat1 - pMat2) / (pS * pS));
            partialS += (pMat1 - pMat2) * ratio;
            partialS_mu += ratio * (pD - pAP);
        }

        float s = partialS.scalar(), s_mu = partialS_mu.scalar();
        for (unsigned j = i.value(); j < size; ++j)
        {
            float ratio = (mat1[j] - mat2[j]) / (S[j] * S[j]);
            s += (mat1[j] - mat2[j]) * ratio;
            s_mu += ratio * (D[j] - AP[j]);
        }
        return AlphaParameters(s, s_mu);
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

// PERFORMANCE_CRITICAL
AlphaParameters DenseGibbsSampler::alphaParametersWithChange(unsigned row,
unsigned col, float ch)
{
    unsigned size = mDMatrix.nRow();
    const float *D = mDMatrix.getCol(row).ptr();
    const float *S = mSMatrix.getCol(row).ptr();
    const float *AP = mAPMatrix.getCol(row).ptr();
    const float *mat = mOtherMatrix->getCol(col).ptr();

    gaps::simd::PackedFloat pCh(ch);
    gaps::simd::PackedFloat pMat, pD, pAP, pS;
    gaps::simd::PackedFloat partialS(0.f), partialS_mu(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        gaps::simd::PackedFloat ratio(pMat / (pS * pS));
        partialS += pMat * ratio;
        partialS_mu += ratio * (pD - (pAP + pCh * pMat));
    }

    float s = partialS.scalar(), s_mu = partialS_mu.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        float ratio = mat[j] / (S[j] * S[j]);
        s += mat[j] * ratio;
        s_mu += ratio * (D[j] - (AP[j] + ch * mat[j]));
    }
    return AlphaParameters(s, s_mu);

}

// PERFORMANCE_CRITICAL
void DenseGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->getCol(col).ptr();
    float *ap = mAPMatrix.getCol(row).ptr();
    unsigned size = mAPMatrix.nRow();

    gaps::simd::PackedFloat pOther, pAP;
    gaps::simd::Index i(0);
    gaps::simd::PackedFloat pDelta(delta);
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

Archive& operator<<(Archive &ar, DenseGibbsSampler &s)
{
    ar << s.mMatrix << s.mDomain << s.mAlpha << s.mLambda << s.mMaxGibbsMass
        << s.mAnnealingTemp << s.mNumPatterns << s.mNumBins << s.mBinLength;
    return ar;
}

Archive& operator>>(Archive &ar, DenseGibbsSampler &s)
{
    ar >> s.mMatrix >> s.mDomain >> s.mAlpha >> s.mLambda >> s.mMaxGibbsMass
        >> s.mAnnealingTemp >> s.mNumPatterns >> s.mNumBins >> s.mBinLength;
    return ar;
}