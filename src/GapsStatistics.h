#ifndef __COGAPS_GAPS_STATISTICS_H__
#define __COGAPS_GAPS_STATISTICS_H__

#include "gibbs_sampler/DenseNormalModel.h"
#include "gibbs_sampler/SparseNormalModel.h"
#include "math/Math.h"
#include "math/MatrixMath.h"
#include "math/VectorMath.h"
#include "data_structures/Matrix.h"
#include "utils/GapsAssert.h"

#define GAPS_SQ(x) ((x) * (x))

class Archive;

class GapsStatistics
{
public:
    GapsStatistics(unsigned nGenes, unsigned nSamples, unsigned nPatterns);
    template <class DataModel>
    void update(const DataModel &AModel, const DataModel &PModel);
    template <class DataModel>
    void updatePump(const DataModel &AModel);
    template <class DataModel>
    void takeSnapshot(GapsAlgorithmPhase whichPhase, const DataModel &AModel, const DataModel &PModel);
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
    float meanChiSq(const DenseNormalModel &model) const;
    float meanChiSq(const SparseNormalModel &model) const;
    const std::vector<Matrix>& getEquilibrationSnapshots(char whichMatrix) const;
    const std::vector<Matrix>& getSamplingSnapshots(char whichMatrix) const;
    friend Archive& operator<<(Archive &ar, const GapsStatistics &stat);
    friend Archive& operator>>(Archive &ar, GapsStatistics &stat);
private:
    Matrix mAMeanMatrix;
    Matrix mAStdMatrix;
    Matrix mPMeanMatrix;
    Matrix mPStdMatrix;
    Matrix mPumpMatrix;
    std::vector<Matrix> mEquilibrationSnapshotsA;
    std::vector<Matrix> mEquilibrationSnapshotsP;
    std::vector<Matrix> mSamplingSnapshotsA;
    std::vector<Matrix> mSamplingSnapshotsP;
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
        statMatrix->operator()(i, maxIndices[i])++;
    }
}

// implemented the same as unique for now - TODO
template <class FactorMatrix>
inline void pumpMatrixCutThreshold(const FactorMatrix &AMatrix, Matrix *statMatrix)
{
    GAPS_ASSERT(statMatrix->nRow() == AMatrix.nRow());
    GAPS_ASSERT(statMatrix->nCol() == AMatrix.nCol());
    std::vector<float> maxValues(AMatrix.nRow(), 0.f);
    std::vector<unsigned> maxIndices(AMatrix.nRow(), 0);
    for (unsigned j = 0; j < AMatrix.nCol(); ++j) // col-major ordering
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

template <class DataModel>
void GapsStatistics::updatePump(const DataModel &AModel)
{
    ++mPumpUpdates;
    if (mPumpThreshold == PUMP_UNIQUE)
    {
        pumpMatrixCutThreshold(AModel.mMatrix, &mPumpMatrix);
    }
    else
    {
        pumpMatrixUniqueThreshold(AModel.mMatrix, &mPumpMatrix);
    }
}

// TODO is there precision loss? use double?
template <class DataModel>
void GapsStatistics::update(const DataModel &AModel, const DataModel &PModel)
{
    GAPS_ASSERT(mNumPatterns == AModel.mMatrix.nCol());
    GAPS_ASSERT(mNumPatterns == PModel.mMatrix.nCol());
    ++mStatUpdates;
    for (unsigned j = 0; j < mNumPatterns; ++j)
    {
        float norm = gaps::max(PModel.mMatrix.getCol(j));
        norm = (norm == 0.f) ? 1.f : norm;
        Vector quot(PModel.mMatrix.getCol(j) / norm);
        mPMeanMatrix.getCol(j) += quot;
        mPStdMatrix.getCol(j) += gaps::elementSq(quot);
        Vector prod(AModel.mMatrix.getCol(j) * norm);
        mAMeanMatrix.getCol(j) += prod;
        mAStdMatrix.getCol(j) += gaps::elementSq(prod);
        GAPS_ASSERT(norm > 0.f);
        GAPS_ASSERT(gaps::min(quot) >= 0.f);
        GAPS_ASSERT(gaps::min(prod) >= 0.f);
    }
}

template <class DataModel>
void GapsStatistics::takeSnapshot(GapsAlgorithmPhase whichPhase, const DataModel &AModel,
const DataModel &PModel)
{
    if (whichPhase == GAPS_EQUILIBRATION_PHASE)
    {
        mEquilibrationSnapshotsA.push_back(AModel.mMatrix.getMatrix());
        mEquilibrationSnapshotsP.push_back(PModel.mMatrix.getMatrix());
    }
    else if (whichPhase == GAPS_SAMPLING_PHASE)
    {
        mSamplingSnapshotsA.push_back(AModel.mMatrix.getMatrix());
        mSamplingSnapshotsP.push_back(PModel.mMatrix.getMatrix());
    }
}

#endif // __COGAPS_GAPS_STATISTICS_H__