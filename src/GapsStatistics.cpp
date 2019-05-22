#include "GapsStatistics.h"
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

float GapsStatistics::meanChiSq(const GibbsSampler<DenseStorage> &PSampler) const
{
    float chisq = 0.f;
    for (unsigned i = 0; i < PSampler.mDMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < PSampler.mDMatrix.nCol(); ++j)
        {
            float m = 0.f;
            for (unsigned k = 0; k < mAMeanMatrix.nCol(); ++k)
            {
                m += mAMeanMatrix(i,k) * mPMeanMatrix(j,k);
            }
            m /= GAPS_SQ(static_cast<float>(mStatUpdates));

            float d = PSampler.mDMatrix(i,j);
            float s = PSampler.mSMatrix(i,j);
            chisq += GAPS_SQ(d - m) / GAPS_SQ(s);
        }
    }
    return chisq;
}

float GapsStatistics::meanChiSq(const GibbsSampler<SparseStorage> &PSampler) const
{
    float chisq = 0.f;
    for (unsigned i = 0; i < PSampler.mDMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < PSampler.mDMatrix.nCol(); ++j)
        {
            float m = 0.f;
            for (unsigned k = 0; k < mAMeanMatrix.nCol(); ++k)
            {
                m += mAMeanMatrix(i,k) * mPMeanMatrix(j,k);
            }
            m /= GAPS_SQ(static_cast<float>(mStatUpdates));

            float d = PSampler.mDMatrix.getCol(j).at(i);
            float s = gaps::max(d * 0.1f, 0.1f);
            chisq += GAPS_SQ(d - m) / GAPS_SQ(s);
        }
    }
    return chisq;
}

float GapsStatistics::meanChiSq(const SingleThreadedGibbsSampler<DenseStorage> &PSampler) const
{
    float chisq = 0.f;
    for (unsigned i = 0; i < PSampler.mDMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < PSampler.mDMatrix.nCol(); ++j)
        {
            float m = 0.f;
            for (unsigned k = 0; k < mAMeanMatrix.nCol(); ++k)
            {
                m += mAMeanMatrix(i,k) * mPMeanMatrix(j,k);
            }
            m /= GAPS_SQ(static_cast<float>(mStatUpdates));

            float d = PSampler.mDMatrix(i,j);
            float s = PSampler.mSMatrix(i,j);
            chisq += GAPS_SQ(d - m) / GAPS_SQ(s);
        }
    }
    return chisq;
}

float GapsStatistics::meanChiSq(const SingleThreadedGibbsSampler<SparseStorage> &PSampler) const
{
    float chisq = 0.f;
    for (unsigned i = 0; i < PSampler.mDMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < PSampler.mDMatrix.nCol(); ++j)
        {
            float m = 0.f;
            for (unsigned k = 0; k < mAMeanMatrix.nCol(); ++k)
            {
                m += mAMeanMatrix(i,k) * mPMeanMatrix(j,k);
            }
            m /= GAPS_SQ(static_cast<float>(mStatUpdates));

            float d = PSampler.mDMatrix.getCol(j).at(i);
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

