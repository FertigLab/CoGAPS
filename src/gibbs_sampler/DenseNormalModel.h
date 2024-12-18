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
public:
    template <class DataType>
    DenseNormalModel(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params, float alpha, float maxGibbsMass);
    template <class DataType>
    void setUncertainty(const DataType &unc, bool transpose, bool subsetRows,
        const GapsParameters &params);
    void setMatrix(const Matrix &mat);
    void setAnnealingTemp(float temp);
    void sync(const DenseNormalModel &model, unsigned nThreads=1);
    void extraInitialization();
    float chiSq() const;
    float dataSparsity() const;
    friend Archive& operator<<(Archive &ar, const DenseNormalModel &m);
    friend Archive& operator>>(Archive &ar, DenseNormalModel &m);
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
    DenseNormalModel(const DenseNormalModel&); // = delete (no c++11)
    DenseNormalModel& operator=(const DenseNormalModel&); // = delete (no c++11)
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);
    void updateAPMatrix(unsigned row, unsigned col, float delta);

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
    std::cout<<" trnsp "<<transpose <<std::endl;
    std::cout << std::scientific<<std::setprecision(40);
    std::cout<<" mean new uncert mtrx "<<gaps::mean(mSMatrix)  <<std::endl;
    std::cout<<" min new uncert mtrx "<<gaps::min(mSMatrix)  <<std::endl;
    std::cout<<" max new uncert mtrx "<<gaps::max(mSMatrix)  <<std::endl;

    mLambda = alpha * std::sqrt(nPatterns() / meanD);
    mMaxGibbsMass = mMaxGibbsMass / mLambda;

    if (gaps::max(mDMatrix) > 50.f)
    {
        gaps_printf("\nWarning: Large values detected, is data log transformed?\n");
    }
    mSMatrix=gaps::pmax(mDMatrix, factor, mLambda); 
    //we suppose that mLambda estimates an atom size expectation
    //so msMatrix[i,j] is max(mDMatrix[i,j]*factor, atom size)
    // it was Tom's mSMatrix.pad(1.f); // so that SIMD operations don't divide by zero
    std::cout<<"mLambda "<<mLambda<<std::endl;
    std::cout<<"min - temp value "<<gaps::min(gaps::pmax(mDMatrix, factor, mLambda))  <<std::endl;
    std::cout<<"After pmax mean uncert mtrx "<<gaps::mean(mSMatrix)  <<std::endl;
    std::cout<<"After pmax min uncert mtrx "<<gaps::min(mSMatrix)  <<std::endl;
    std::cout<<"After pmax max uncert mtrx "<<gaps::max(mSMatrix)  <<std::endl;
}

template <class DataType>
void DenseNormalModel::setUncertainty(const DataType &unc, bool transpose,
bool subsetRows, const GapsParameters &params)
{
    mSMatrix = Matrix(unc, transpose, subsetRows, params.dataIndicesSubset);
    mSMatrix.pad(1.f); // so that SIMD operations don't divide by zero
}

#endif // __COGAPS_DENSE_STORAGE_POLICY_H__
