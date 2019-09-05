#include "GapsStatistics.h"
#include "utils/Archive.h"
#include "math/Math.h"

GapsStatistics::GapsStatistics(unsigned nGenes, unsigned nSamples, unsigned nPatterns)
    :
mAMeanMatrix(nGenes, nPatterns), mAStdMatrix(nGenes, nPatterns),
mPMeanMatrix(nSamples, nPatterns), mPStdMatrix(nSamples, nPatterns),
mPumpMatrix(nGenes, nPatterns), mPumpThreshold(PUMP_CUT), mStatUpdates(0),
mNumPatterns(nPatterns), mPumpUpdates(0)
{}

Matrix GapsStatistics::Amean() const
{
#ifdef GAPS_DEBUG
    gaps_printf("max value of Amean: %f\n", gaps::max(mAMeanMatrix));
    gaps_printf("min value of Amean: %f\n", gaps::min(mAMeanMatrix));
#endif
    GAPS_ASSERT(gaps::min(mAMeanMatrix / static_cast<float>(mStatUpdates)) >= 0.f);
    return mAMeanMatrix / static_cast<float>(mStatUpdates);
}

Matrix GapsStatistics::Asd() const
{
    Matrix mat(mAStdMatrix.nRow(), mAStdMatrix.nCol());
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            float meanTerm = GAPS_SQ(mAMeanMatrix(i,j)) / static_cast<float>(mStatUpdates);
            float numer = gaps::max(0.f, mAStdMatrix(i,j) - meanTerm);
            mat(i,j) = std::sqrt(numer / (static_cast<float>(mStatUpdates) - 1.f));
        }
    }
    return mat;
}

Matrix GapsStatistics::Pmean() const
{
#ifdef GAPS_DEBUG
    gaps_printf("max value of Pmean: %f\n", gaps::max(mPMeanMatrix));
    gaps_printf("min value of Pmean: %f\n", gaps::min(mPMeanMatrix));
#endif
    GAPS_ASSERT(gaps::min(mPMeanMatrix / static_cast<float>(mStatUpdates)) >= 0.f);
    return mPMeanMatrix / static_cast<float>(mStatUpdates);
}

Matrix GapsStatistics::Psd() const
{
    Matrix mat(mPStdMatrix.nRow(), mPStdMatrix.nCol());
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            float meanTerm = GAPS_SQ(mPMeanMatrix(i,j)) / static_cast<float>(mStatUpdates);
            float numer = gaps::max(0.f, mPStdMatrix(i,j) - meanTerm);
            mat(i,j) = std::sqrt(numer / (static_cast<float>(mStatUpdates) - 1.f));
        }
    }
    return mat;
}

float GapsStatistics::meanChiSq(const DenseNormalModel &model) const
{
    GAPS_ASSERT(model.mDMatrix.nRow() == mAMeanMatrix.nRow());
    GAPS_ASSERT(model.mDMatrix.nCol() == mPMeanMatrix.nRow());

    float chisq = 0.f;
    for (unsigned i = 0; i < model.mDMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < model.mDMatrix.nCol(); ++j)
        {
            float m = 0.f;
            for (unsigned k = 0; k < mAMeanMatrix.nCol(); ++k)
            {
                m += mAMeanMatrix(i,k) * mPMeanMatrix(j,k);
            }
            m /= GAPS_SQ(static_cast<float>(mStatUpdates));

            float d = model.mDMatrix(i,j);
            float s = model.mSMatrix(i,j);
            chisq += GAPS_SQ(d - m) / GAPS_SQ(s);
        }
    }
    return chisq;
}

float GapsStatistics::meanChiSq(const SparseNormalModel &model) const
{
    GAPS_ASSERT(model.mDMatrix.nRow() == mAMeanMatrix.nRow());
    GAPS_ASSERT(model.mDMatrix.nCol() == mPMeanMatrix.nRow());

    float chisq = 0.f;
    for (unsigned i = 0; i < model.mDMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < model.mDMatrix.nCol(); ++j)
        {
            float m = 0.f;
            for (unsigned k = 0; k < mAMeanMatrix.nCol(); ++k)
            {
                m += mAMeanMatrix(i,k) * mPMeanMatrix(j,k);
            }
            m /= GAPS_SQ(static_cast<float>(mStatUpdates));

            float d = model.mDMatrix.getCol(j).at(i);
            float s = gaps::max(d * 0.1f, 0.1f);
            chisq += GAPS_SQ(d - m) / GAPS_SQ(s);
        }
    }
    return chisq;
}

Matrix GapsStatistics::pumpMatrix() const
{
    float denom = mPumpUpdates != 0 ? static_cast<float>(mPumpUpdates) : 1.f;
    return mPumpMatrix / denom;
}

Matrix GapsStatistics::meanPattern() const
{
    Matrix mat(mAMeanMatrix.nRow(), mAMeanMatrix.nCol());
    if (mPumpThreshold == PUMP_UNIQUE)
    {
        pumpMatrixCutThreshold(Amean(), &mat);
    }
    else
    {
        pumpMatrixUniqueThreshold(Amean(), &mat);
    }   
    return mat;
}

void GapsStatistics::addChiSq(float chisq)
{
    mChisqHistory.push_back(chisq);
}

void GapsStatistics::addAtomCount(unsigned atomA, unsigned atomP)
{
    mAtomHistoryA.push_back(atomA);
    mAtomHistoryP.push_back(atomP);
}

std::vector<float> GapsStatistics::chisqHistory() const
{
    return mChisqHistory;
}

std::vector<unsigned> GapsStatistics::atomHistory(char m) const
{
    return m == 'A' ? mAtomHistoryA : mAtomHistoryP;
}

const std::vector<Matrix>& GapsStatistics::getEquilibrationSnapshots(char whichMatrix) const
{
    return (whichMatrix == 'A') ? mEquilibrationSnapshotsA : mEquilibrationSnapshotsP;
}

const std::vector<Matrix>& GapsStatistics::getSamplingSnapshots(char whichMatrix) const
{
    return (whichMatrix == 'A') ? mSamplingSnapshotsA : mSamplingSnapshotsP;
}

Archive& operator<<(Archive &ar, const GapsStatistics &stat)
{
    ar << stat.mAMeanMatrix << stat.mAStdMatrix << stat.mPMeanMatrix
        << stat.mPStdMatrix << stat.mStatUpdates << stat.mNumPatterns;
    return ar;
}

Archive& operator>>(Archive &ar, GapsStatistics &stat)
{
    ar >> stat.mAMeanMatrix >> stat.mAStdMatrix >> stat.mPMeanMatrix
        >> stat.mPStdMatrix >> stat.mStatUpdates >> stat.mNumPatterns;
    return ar;
}

