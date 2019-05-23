#include "DenseNormalWithUncertaintyModel.h"
#include "../math/Math.h"

#define GAPS_SQ(x) ((x) * (x))

void DenseNormalWithUncertaintyModel::setMatrix(const Matrix &mat)
{
    mMatrix = mat;
}

void DenseNormalWithUncertaintyModel::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

// copy transpose of other AP matrix
void DenseNormalWithUncertaintyModel::sync(const DenseNormalWithUncertaintyModel &model, unsigned nThreads)
{
    GAPS_ASSERT(model.mAPMatrix.nRow() == mAPMatrix.nCol());
    GAPS_ASSERT(model.mAPMatrix.nCol() == mAPMatrix.nRow());

    unsigned nc = model.mAPMatrix.nCol();
    unsigned nr = model.mAPMatrix.nRow();

    #pragma omp parallel for num_threads(nThreads)
    for (unsigned j = 0; j < nc; ++j)
    {
        for (unsigned i = 0; i < nr; ++i)
        {
            mAPMatrix(j,i) = model.mAPMatrix(i,j);
        }
    }
    mOtherMatrix = &(model.mMatrix); // update pointer
    GAPS_ASSERT(mOtherMatrix->nCol() == mMatrix.nCol());
}

void DenseNormalWithUncertaintyModel::extraInitialization()
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

float DenseNormalWithUncertaintyModel::chiSq() const
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

float DenseNormalWithUncertaintyModel::dataSparsity() const
{
    return gaps::sparsity(mDMatrix);
}

uint64_t DenseNormalWithUncertaintyModel::nElements() const
{
    return mDMatrix.nRow() * mMatrix.nCol();
}

uint64_t DenseNormalWithUncertaintyModel::nPatterns() const
{
    return mMatrix.nCol();
}

float DenseNormalWithUncertaintyModel::lambda() const
{
    return mLambda;
}

bool DenseNormalWithUncertaintyModel::canUseGibbs(unsigned col) const
{
    return !gaps::isVectorZero(mOtherMatrix->getCol(col));
}

bool DenseNormalWithUncertaintyModel::canUseGibbs(unsigned c1, unsigned c2) const
{
    return canUseGibbs(c1) || canUseGibbs(c2);
}

void DenseNormalWithUncertaintyModel::changeMatrix(unsigned row, unsigned col, float delta)
{
    mMatrix(row, col) += delta;
    updateAPMatrix(row, col, delta);

    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
}

void DenseNormalWithUncertaintyModel::safelyChangeMatrix(unsigned row, unsigned col, float delta)
{
    float newVal = gaps::max(mMatrix(row, col) + delta, 0.f);
    updateAPMatrix(row, col, newVal - mMatrix(row, col));
    mMatrix(row, col) = newVal;

    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
}

float DenseNormalWithUncertaintyModel::deltaLogLikelihood(unsigned r1, unsigned c1, unsigned r2,
unsigned c2, float mass)
{
    AlphaParameters alpha = alphaParameters(r1, c1, r2, c2) * mAnnealingTemp;
    return -1.f * mass * (alpha.s_mu + alpha.s * mass / 2.f);
}

OptionalFloat DenseNormalWithUncertaintyModel::sampleBirth(unsigned row, unsigned col, GapsRng *rng)
{
    AlphaParameters alpha = alphaParameters(row, col) * mAnnealingTemp;
    return gibbsMass(alpha, 0.f, mMaxGibbsMass, rng, mLambda);
}

OptionalFloat DenseNormalWithUncertaintyModel::sampleDeathAndRebirth(unsigned row, unsigned col,
float delta, GapsRng *rng)
{
    AlphaParameters alpha = alphaParametersWithChange(row, col, delta) * mAnnealingTemp;
    OptionalFloat mass = gibbsMass(alpha, 0.f, mMaxGibbsMass, rng, mLambda);
    if (mass.hasValue())
    {
        float deltaLL = mass.value() * (alpha.s_mu - alpha.s * mass.value() / 2.f);
        if (std::log(rng->uniform()) < deltaLL)
        {
            return mass;
        }
    }
    return OptionalFloat();
}

OptionalFloat DenseNormalWithUncertaintyModel::sampleExchange(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2, GapsRng *rng)
{
    AlphaParameters alpha = alphaParameters(r1, c1, r2, c2) * mAnnealingTemp;
    return gibbsMass(alpha, -m1, m2, rng);
}

// PERFORMANCE_CRITICAL
AlphaParameters DenseNormalWithUncertaintyModel::alphaParameters(unsigned row, unsigned col)
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
    return AlphaParameters(partialS.scalar(), partialS_mu.scalar());
}

// PERFORMANCE_CRITICAL
AlphaParameters DenseNormalWithUncertaintyModel::alphaParameters(unsigned r1, unsigned c1,
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
AlphaParameters DenseNormalWithUncertaintyModel::alphaParametersWithChange(unsigned row,
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
void DenseNormalWithUncertaintyModel::updateAPMatrix(unsigned row, unsigned col, float delta)
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

Archive& operator<<(Archive &ar, const DenseNormalWithUncertaintyModel &mod)
{
    ar << mod.mMatrix;
    return ar;
}

Archive& operator>>(Archive &ar, DenseNormalWithUncertaintyModel &mod)
{
    ar >> mod.mMatrix;
    return ar;
}

