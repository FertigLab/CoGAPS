#include "GapsStatistics.h"
#include "math/Math.h"

GapsStatistics::GapsStatistics(unsigned nGenes, unsigned nSamples, unsigned nPatterns)
    :
mAMeanMatrix(nGenes, nPatterns), mAStdMatrix(nGenes, nPatterns),
mPMeanMatrix(nSamples, nPatterns), mPStdMatrix(nSamples, nPatterns),
mStatUpdates(0), mNumPatterns(nPatterns)
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

float GapsStatistics::meanChiSq(const DenseGibbsSampler &PSampler) const
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
            chisq += GAPS_SQ((PSampler.mDMatrix(i,j) - m) / PSampler.mSMatrix(i,j));
        }
    }
    return chisq;
}

float GapsStatistics::meanChiSq(const SparseGibbsSampler &PSampler) const
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
            chisq += GAPS_SQ((d - m) / s);
        }
    }
    return chisq;
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

