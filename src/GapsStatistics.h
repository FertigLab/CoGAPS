#ifndef __COGAPS_GAPS_STATISTICS_H__
#define __COGAPS_GAPS_STATISTICS_H__

#include "GapsParameters.h"
#include "gibbs_sampler/GibbsSampler.h"
#include "gibbs_sampler/DenseStoragePolicy.h"
#include "gibbs_sampler/SparseStoragePolicy.h"
#include "math/Math.h"
#include "math/MatrixMath.h"
#include "math/VectorMath.h"
#include "data_structures/Matrix.h"

#define GAPS_SQ(x) ((x) * (x))

class GapsStatistics
{
public:

    GapsStatistics(unsigned nGenes, unsigned nSamples, unsigned nPatterns);

    template <class Sampler>
    void update(const Sampler &ASampler, const Sampler &PSampler);

    template <class Sampler>
    void updatePump(const Sampler &ASampler);

    Matrix Amean() const;
    Matrix Pmean() const;
    Matrix Asd() const;
    Matrix Psd() const;

    Matrix pumpMatrix() const;
    Matrix meanPattern() const;

    void addChiSq(float chisq);
    void addAtomCount(unsigned atomA, unsigned atomP);

    std::vector<float> chisqHistory() const;
    std::vector<unsigned> atomHistory(char m) const;
    
    float meanChiSq(const GibbsSampler<DenseStorage> &PSampler) const;
    float meanChiSq(const GibbsSampler<SparseStorage> &PSampler) const;

    // serialization
    friend Archive& operator<<(Archive &ar, const GapsStatistics &stat);
    friend Archive& operator>>(Archive &ar, GapsStatistics &stat);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    Matrix mAMeanMatrix;
    Matrix mAStdMatrix;
    Matrix mPMeanMatrix;
    Matrix mPStdMatrix;
    Matrix mPumpMatrix;
   
    std::vector<float> mChisqHistory;
    std::vector<unsigned> mAtomHistoryA;
    std::vector<unsigned> mAtomHistoryP;

    PumpThreshold mPumpThreshold;

    unsigned mStatUpdates;
    unsigned mNumPatterns;
    unsigned mPumpUpdates;
};

template <class FactorMatrix>
inline void pumpMatrixUniqueThreshold(const FactorMatrix &AMatrix, Matrix *statMatrix)
{
    GAPS_ASSERT(statMatrix->nRow() == AMatrix.nRow());
    GAPS_ASSERT(statMatrix->nCol() == AMatrix.nCol());

    std::vector<float> maxValues(AMatrix.nRow(), 0.f);
    std::vector<unsigned> maxIndices(AMatrix.nRow(), 0);
    for (unsigned j = 0; j < AMatrix.nCol(); ++j)
    {
        for (unsigned i = 0; i < AMatrix.nRow(); ++i)
        {
            if (maxValues[i] < AMatrix(i,j))
            {
                maxValues[i] = AMatrix(i,j);
                maxIndices[i] = j;
            }
        }
    }

    for (unsigned i = 0; i < AMatrix.nRow(); ++i)
    {
        statMatrix->operator()(i,maxIndices[i])++;
    }
}

// implemented the same as unique for now - TODO
template <class FactorMatrix>
inline void pumpMatrixCutThreshold(const FactorMatrix &AMatrix, Matrix *statMatrix)
{
    GAPS_ASSERT(statMatrix->nRow() == AMatrix.nRow());
    GAPS_ASSERT(statMatrix->nCol() == AMatrix.nCol());

    // we need to access data in columns due to matrix data layout
    std::vector<float> maxValues(AMatrix.nRow(), 0.f);
    std::vector<unsigned> maxIndices(AMatrix.nRow(), 0);
    for (unsigned j = 0; j < AMatrix.nCol(); ++j)
    {
        for (unsigned i = 0; i < AMatrix.nRow(); ++i)
        {
            if (maxValues[i] < AMatrix(i,j))
            {
                maxValues[i] = AMatrix(i,j);
                maxIndices[i] = j;
            }
        }
    }

    for (unsigned i = 0; i < AMatrix.nRow(); ++i)
    {
        statMatrix->operator()(i,maxIndices[i])++;
    }
}

template <class Sampler>
void GapsStatistics::updatePump(const Sampler &ASampler)
{
    ++mPumpUpdates;
    if (mPumpThreshold == PUMP_UNIQUE)
    {
        pumpMatrixCutThreshold(ASampler.mMatrix, &mPumpMatrix);
    }
    else
    {
        pumpMatrixUniqueThreshold(ASampler.mMatrix, &mPumpMatrix);
    }
}

template <class Sampler>
void GapsStatistics::update(const Sampler &ASampler, const Sampler &PSampler)
{
    mStatUpdates++;

    GAPS_ASSERT(mNumPatterns == ASampler.mMatrix.nCol());
    GAPS_ASSERT(mNumPatterns == PSampler.mMatrix.nCol());

    // precision loss? use double?
    for (unsigned j = 0; j < mNumPatterns; ++j)
    {
        float norm = gaps::max(PSampler.mMatrix.getCol(j));
        norm = norm == 0.f ? 1.f : norm;
        GAPS_ASSERT(norm > 0.f);

        Vector quot(PSampler.mMatrix.getCol(j) / norm);
        GAPS_ASSERT(gaps::min(quot) >= 0.f);
        mPMeanMatrix.getCol(j) += quot;
        mPStdMatrix.getCol(j) += gaps::elementSq(quot);

        Vector prod(ASampler.mMatrix.getCol(j) * norm);
        GAPS_ASSERT(gaps::min(prod) >= 0.f);
        mAMeanMatrix.getCol(j) += prod;
        mAStdMatrix.getCol(j) += gaps::elementSq(prod);
    }
}

#endif