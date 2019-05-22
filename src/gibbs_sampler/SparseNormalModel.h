#ifndef __COGAPS_SPARSE_NORMAL_MODEL_H__
#define __COGAPS_SPARSE_NORMAL_MODEL_H__

#include "AlphaParameters.h"
#include "../GapsParameters.h"
#include "../data_structures/HybridMatrix.h"
#include "../data_structures/SparseMatrix.h"
#include "../math/MatrixMath.h"

class GapsStatistics;
class Archive;

class SparseNormalModel
{
public:

    friend class GapsStatistics;

    template <class DataType>
    SparseNormalModel(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params, float alpha);

    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params);

    void setMatrix(const Matrix &mat);
    void setAnnealingTemp(float temp);

    void sync(const SparseNormalModel &sampler, unsigned nThreads=1);
    void extraInitialization();

    float chiSq() const;
    float dataSparsity() const;

    friend Archive& operator<<(Archive &ar, const SparseNormalModel &s);
    friend Archive& operator>>(Archive &ar, SparseNormalModel &s);

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

private:

    void generateLookupTables();

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);

    SparseNormalModel(const SparseNormalModel&); // = delete (no c++11)
    SparseNormalModel& operator=(const SparseNormalModel&); // = delete (no c++11)

    SparseMatrix mDMatrix; // samples by genes for A, genes by samples for P
    HybridMatrix mMatrix; // genes by patterns for A, samples by patterns for P
    const HybridMatrix *mOtherMatrix; // pointer to P if this is A, and vice versa

    Matrix mZ2;
    Vector mZ1;

    float mBeta;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    float mLambda;
};

template <class DataType>
SparseNormalModel::SparseNormalModel(const DataType &data, bool transpose,
bool subsetRows, const GapsParameters &params, float alpha)
    :
mDMatrix(data, transpose, subsetRows, params.dataIndicesSubset),
mMatrix(mDMatrix.nCol(), params.nPatterns),
mOtherMatrix(NULL),
mZ2(params.nPatterns, params.nPatterns),
mZ1(params.nPatterns),
mBeta(100.f),
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
}

// required function for the GibbsSampler interface
template <class DataType>
void SparseNormalModel::setUncertainty(const DataType &data, bool transpose,
bool subsetRows, const GapsParameters &params)
{
    // nop - SparseGibbsSampler assumes default uncertainty always
}

#endif // __COGAPS_SPARSE_NORMAL_MODEL_H__