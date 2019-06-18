#ifndef __COGAPS_DENSE_NORMAL_MODEL_H__
#define __COGAPS_DENSE_NORMAL_MODEL_H__

#include "AlphaParameters.h"
#include "../GapsParameters.h"
#include "../data_structures/Matrix.h"
#include "../math/MatrixMath.h"

class GapsStatistics;
class Archive;

class DenseNormalModel
{
public:

    friend class GapsStatistics;

    template <class DataType>
    DenseNormalModel(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params, float alpha);

    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params);

    void setMatrix(const Matrix &mat);
    void setAnnealingTemp(float temp);

    void sync(const DenseNormalModel &model, unsigned nThreads=1);
    void extraInitialization();

    float chiSq() const;
    float dataSparsity() const;

    friend Archive& operator<<(Archive &ar, const DenseNormalModel &s);
    friend Archive& operator>>(Archive &ar, DenseNormalModel &s);

protected:

    uint64_t nElements() const;
    uint64_t nPatterns() const;

    float lambda() const;

    bool canUseGibbs(unsigned col) const;
    bool canUseGibbs(unsigned c1, unsigned c2) const;

    void changeMatrix(unsigned row, unsigned col, float delta);
    void safelyChangeMatrix(unsigned row, unsigned col, float delta);

    float deltaLogLikelihood(unsigned r1, unsigned c1, unsigned r2, unsigned c2, float mass);
    OptionalFloat sampleBirth(unsigned row, unsigned col, GapsRng *rng);
    OptionalFloat sampleDeathAndRebirth(unsigned row, unsigned col, float delta, GapsRng *rng);
    OptionalFloat sampleExchange(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2, GapsRng *rng);

//private:

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    DenseNormalModel(const DenseNormalModel&); // = delete (no c++11)
    DenseNormalModel& operator=(const DenseNormalModel&); // = delete (no c++11)

    Matrix mDMatrix; // samples by genes for A, genes by samples for P
    Matrix mMatrix; // genes by patterns for A, samples by patterns for P
    const Matrix *mOtherMatrix; // pointer to P if this is A, and vice versa

    Matrix mSMatrix; // uncertainty values for each data point
    Matrix mAPMatrix; // cached product of A and P

    float mMaxGibbsMass;
    float mAnnealingTemp;
    float mLambda;
};

template <class DataType>
DenseNormalModel::DenseNormalModel(const DataType &data, bool transpose,
bool subsetRows, const GapsParameters &params, float alpha)
    :
mDMatrix(data, transpose, subsetRows, params.dataIndicesSubset),
mMatrix(mDMatrix.nCol(), params.nPatterns),
mOtherMatrix(NULL),
mSMatrix(gaps::pmax(mDMatrix, 0.1f)),
mAPMatrix(mDMatrix.nRow(), mDMatrix.nCol()),
mMaxGibbsMass(100.f),
mAnnealingTemp(1.f),
mLambda(0.f)
{
    float meanD = params.singleCell ? gaps::nonZeroMean(mDMatrix) : gaps::mean(mDMatrix);
    mLambda = alpha * std::sqrt(nPatterns() / meanD);
    mMaxGibbsMass = mMaxGibbsMass / mLambda;

    if (gaps::max(mDMatrix) > 50.f)
    {
        gaps_printf("\nWarning: Large values detected, is data log transformed?\n");
    }
    mSMatrix.pad(1.f); // so that SIMD operations don't divide by zero
}

template <class DataType>
void DenseNormalModel::setUncertainty(const DataType &unc, bool transpose,
bool subsetRows, const GapsParameters &params)
{
    mSMatrix = Matrix(unc, transpose, subsetRows, params.dataIndicesSubset);
    mSMatrix.pad(1.f); // so that SIMD operations don't divide by zero
}

#endif // __COGAPS_DENSE_STORAGE_POLICY_H__