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
    GAPS_ASSERT_MSG(mOtherMatrix->nCol() == mMatrix.nCol(),
        mOtherMatrix->nCol() << " != " << mMatrix.nCol());
}

void DenseGibbsSampler::extraInitialization()
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

    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
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
    for (gaps::simd::Index i(0); i < size; ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        gaps::simd::PackedFloat ratio(pMat / (pS * pS));
        partialS += pMat * ratio;
        partialS_mu += ratio * (pD - pAP);
    }
    
    float s = partialS.scalar();
    float s_mu = partialS_mu.scalar();

#if 0 // optional debug section, set to 1 to enable, 0 to disable
#ifdef GAPS_DEBUG
    float debug_s = 0.f, debug_smu = 0.f;
    for (unsigned j = 0; j < size; ++j)
    {
        float ratio = mat[j] / (S[j] * S[j]);
        debug_s += mat[j] * ratio;
        debug_smu += ratio * (D[j] - AP[j]);
    }
   
    float s_denom = debug_s < 10.f ? 10.f : debug_s;
    float smu_denom = std::abs(debug_smu) < 10.f ? 10.f : std::abs(debug_smu);

    const float tolerance = 0.01f;
    GAPS_ASSERT_MSG(std::abs(s - debug_s) / s_denom < tolerance,
        s << " != " << debug_s);
    GAPS_ASSERT_MSG(std::abs(s_mu - debug_smu) / smu_denom < tolerance,
        s_mu << " != " << debug_smu);
#endif
#endif

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
        gaps::simd::PackedFloat packedS(0.f), packedS_mu(0.f);
        for (gaps::simd::Index i(0); i < size; ++i)
        {   
            pMat1.load(mat1 + i);
            pMat2.load(mat2 + i);
            pD.load(D + i);
            pAP.load(AP + i);
            pS.load(S + i);
            gaps::simd::PackedFloat ratio((pMat1 - pMat2) / (pS * pS));
            packedS += (pMat1 - pMat2) * ratio;
            packedS_mu += ratio * (pD - pAP);
        }
        return AlphaParameters(packedS.scalar(), packedS_mu.scalar());
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
    gaps::simd::PackedFloat packedS(0.f), packedS_mu(0.f);
    for (gaps::simd::Index i(0); i < size; ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        gaps::simd::PackedFloat ratio(pMat / (pS * pS));
        packedS += pMat * ratio;
        packedS_mu += ratio * (pD - (pAP + pCh * pMat));
    }
    return AlphaParameters(packedS.scalar(), packedS_mu.scalar());
}

// PERFORMANCE_CRITICAL
void DenseGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->getCol(col).ptr();
    float *ap = mAPMatrix.getCol(row).ptr();
    unsigned size = mAPMatrix.nRow();

    gaps::simd::PackedFloat pOther, pAP;
    gaps::simd::PackedFloat pDelta(delta);
    for (gaps::simd::Index i(0); i < size; ++i)
    {
        pOther.load(other + i);
        pAP.load(ap + i);
        pAP += pDelta * pOther;
        pAP.store(ap + i);
    }
}

Archive& operator<<(Archive &ar, const DenseGibbsSampler &s)
{
    ar << s.mMatrix << s.mDomain << s.mQueue << s.mAlpha << s.mLambda
        << s.mMaxGibbsMass << s.mAnnealingTemp << s.mNumPatterns << s.mNumBins
        << s.mBinLength;
    return ar;
}

Archive& operator>>(Archive &ar, DenseGibbsSampler &s)
{
    ar >> s.mMatrix >> s.mDomain >> s.mQueue >> s.mAlpha >> s.mLambda
        >> s.mMaxGibbsMass >> s.mAnnealingTemp >> s.mNumPatterns >> s.mNumBins
        >> s.mBinLength;
    return ar;
}

#ifdef GAPS_DEBUG
bool DenseGibbsSampler::internallyConsistent() const
{
    for (unsigned j = 0; j < mMatrix.nCol(); ++j)
    {
        for (unsigned i = 0; i < mMatrix.nRow(); ++i)
        {
            if (mMatrix(i,j) < 0.f)
            {
                return false;
            }
        }
    }
    return true;
}
#endif