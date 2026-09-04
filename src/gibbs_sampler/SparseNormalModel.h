#ifndef __COGAPS_SPARSE_NORMAL_MODEL_H__
#define __COGAPS_SPARSE_NORMAL_MODEL_H__

#include "AlphaParameters.h"
#include "../GapsParameters.h"
#include "../data_structures/HybridMatrix.h"
#include "../data_structures/SparseMatrix.h"
#include "../math/MatrixMath.h"
#include "../utils/GapsPrint.h"

#include <cmath>

class GapsStatistics;
class Archive;

class SparseNormalModel
{
public:
    template <class DataType>
    SparseNormalModel(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params, float alpha, float maxGibbsMass);
    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params);
    void setMatrix(const Matrix &mat);
    void setAnnealingTemp(float temp);
    void sync(const SparseNormalModel &model, unsigned nThreads=1);
    void extraInitialization();
    float chiSq() const;
    float dataSparsity() const;
    friend Archive& operator<<(Archive &ar, const SparseNormalModel &m);
    friend Archive& operator>>(Archive &ar, SparseNormalModel &m);
protected:
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
    friend class GapsStatistics;
    SparseNormalModel(const SparseNormalModel&); // = delete (no c++11)
    SparseNormalModel& operator=(const SparseNormalModel&); // = delete (no c++11)
    void generateLookupTables();
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);

    SparseMatrix mDMatrix; // Data Matrix D, samp x genes or genes x samp
    HybridMatrix mMatrix; // A (left mult) or P (right mult) based on mDMatrix
    const HybridMatrix *mOtherMatrix; // pointer to vis-a-vis of mMartix
    Matrix mZ2;
    Vector mZ1;
    float mBeta;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    float mLambda;
};

template <class DataType>
SparseNormalModel::SparseNormalModel(const DataType &data, bool transpose,
bool subsetRows, const GapsParameters &params, float alpha, float maxGibbsMass)
    :
mDMatrix(data, transpose, subsetRows, params.dataIndicesSubset),
mMatrix(mDMatrix.nCol(), params.nPatterns),
mOtherMatrix(NULL),
mZ2(params.nPatterns, params.nPatterns),
mZ1(params.nPatterns),
mBeta(100.f),
mMaxGibbsMass(maxGibbsMass),
mAnnealingTemp(1.f),
mLambda(0.f)
{
    float meanD = gaps::nonZeroMean(mDMatrix);
    mLambda = alpha * std::sqrt(nPatterns() / meanD);
    mMaxGibbsMass = mMaxGibbsMass / mLambda;

    if (gaps::max(mDMatrix) > 50.f)
    {
        gaps_printf("\nWarning: Large values detected, is data log transformed?\n");
    }
}

// required function for the GibbsSampler interface
template <class DataType>
void SparseNormalModel::setUncertainty(const DataType &data, bool transpose, // NOLINT
bool subsetRows, const GapsParameters &params) // NOLINT
{
    // nop - SparseGibbsSampler assumes default uncertainty always
}

#endif // __COGAPS_SPARSE_NORMAL_MODEL_H__