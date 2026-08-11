#ifndef __COGAPS_DENSE_NORMAL_MODEL_H__
#define __COGAPS_DENSE_NORMAL_MODEL_H__

#include "AlphaParameters.h"
#include "../GapsParameters.h"
#include "../data_structures/Matrix.h"
#include "../math/MatrixMath.h"
#include "../utils/GapsPrint.h"

#include <cmath>

class GapsStatistics;
class Archive;

class DenseNormalModel
{
protected:
    friend class GapsStatistics;
    uint64_t nElements() const;
    uint64_t nPatterns() const;
    float annealingTemp() const;
    float lambda() const;
    float maxGibbsMass() const;
    bool canUseGibbs(unsigned col) const;
    bool canUseGibbs(unsigned c1, unsigned c2) const;
    void changeMatrix(unsigned row, unsigned col, float delta);
    void safelyChangeMatrix(unsigned row, unsigned col, float delta);
    float deltaLogLikelihood(unsigned r1, unsigned c1, unsigned r2, unsigned c2, float mass);
    OptionalFloat sampleBirth(unsigned row, unsigned col, GapsRng *rng);
    OptionalFloat sampleDeathAndRebirth(unsigned row, unsigned col, float delta, GapsRng *rng);
    OptionalFloat sampleExchange(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2, GapsRng *rng);
//private: // TODO
    // P means transpose = TRUE
    Matrix mDMatrix; // samples by genes for A, genes by samples for P
    Matrix mMatrix; // genes by patterns for A, samples by patterns for P
    const Matrix *mOtherMatrix; // pointer to P if this is A, and vice versa
    Matrix mSMatrix; // uncertainty values for each data point
    Matrix mAPMatrix; // cached product of A and P
    //GAPS_ASSERT(mMatrix.nRow() == mAPMatrix.nCol());
    
    DenseNormalModel(const DenseNormalModel&); // = delete (no c++11)
    DenseNormalModel& operator=(const DenseNormalModel&); // = delete (no c++11)
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    float mMaxGibbsMass;
    float mAnnealingTemp;
    float mLambda;
public:
    template <class DataType>
    DenseNormalModel(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params, float alpha, float maxGibbsMass);
    template <class DataType>
    void setUncertainty(const DataType &unc, bool transpose, bool subsetRows,
        const GapsParameters &params);
    void setMatrix(const Matrix &mat);
    void setAnnealingTemp(float temp);
    void sync(const DenseNormalModel &model);
    void extraInitialization();
    float chiSq() const;
    float dataSparsity() const;
    const Matrix & APMatrix () const 
    {
        return mAPMatrix;
    };
    const Matrix & MyMatrix () const 
    {
        return mMatrix;
    };
    const Matrix & UMatrix () const 
    {
        return mSMatrix;
    };
    friend Archive& operator<<(Archive &ar, const DenseNormalModel &m);
    friend Archive& operator>>(Archive &ar, DenseNormalModel &m);
};


template <class DataType>
DenseNormalModel::DenseNormalModel(const DataType &data, bool transpose,
bool subsetRows, const GapsParameters &params, float alpha, float maxGibbsMass)
    :
mDMatrix(data, transpose, subsetRows, params.dataIndicesSubset),
mMatrix(mDMatrix.nCol(), params.nPatterns),
mOtherMatrix(NULL),
mSMatrix(mDMatrix.nRow(), mDMatrix.nCol()),
mAPMatrix(mDMatrix.nRow(), mDMatrix.nCol()),
mMaxGibbsMass(maxGibbsMass),
mAnnealingTemp(1.f),
mLambda(0.f)
{
    float meanD = gaps::nonZeroMean(mDMatrix);
    float factor=0.1f; //it is like 42 but for variance

    mLambda = alpha * std::sqrt(nPatterns() / meanD);
    mMaxGibbsMass = mMaxGibbsMass / mLambda;

    if (gaps::max(mDMatrix) > 50.f)
    {
        gaps_printf("\nWarning: Large values detected, is data log transformed?\n");
    }
    // uncertainty model: relative error S = factor*D, floored at factor so that
    // zeros (and any D < 1) get S = factor. This matches the SparseNormalModel
    // assumption (S = D for observed / 1 for zero, scaled by mBeta = 1/factor^2),
    // so sparseOptimization gives the same result as the dense sampler.
    // (mLambda is the atom-size scale used for mMaxGibbsMass only, NOT for S.)
    mSMatrix = gaps::pmax(mDMatrix, factor, factor); // = max(factor*D, factor)
}

template <class DataType>
void DenseNormalModel::setUncertainty(const DataType &unc, bool transpose,
bool subsetRows, const GapsParameters &params)
{
    mSMatrix = Matrix(unc, transpose, subsetRows, params.dataIndicesSubset);
    // Only the SIMD padding may be overwritten -- pad() would set *every* element
    // to 1.f and so discard the uncertainty the caller passed in. Padding lanes
    // get 1.f so that the SIMD loops divide by 1, not by 0.
    mSMatrix.padSIMD(1.f);
}


#endif // __COGAPS_DENSE_STORAGE_POLICY_H__
