#include "GapsStatistics.h"
#include "math/Math.h"

GapsStatistics::GapsStatistics(unsigned nRow, unsigned nCol, unsigned nPatterns)
    :
mAMeanMatrix(nRow, nPatterns), mAStdMatrix(nRow, nPatterns),
mPMeanMatrix(nCol, nPatterns), mPStdMatrix(nCol, nPatterns),
mStatUpdates(0), mNumPatterns(nPatterns)
{}

Matrix GapsStatistics::Amean() const
{
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
    Matrix A(mAMeanMatrix / static_cast<float>(mStatUpdates));
    Matrix P(mPMeanMatrix / static_cast<float>(mStatUpdates));
    
    float chisq = 0.f;
    for (unsigned i = 0; i < PSampler.mDMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < PSampler.mDMatrix.nCol(); ++j)
        {
            float m = 0.f;
            for (unsigned k = 0; k < A.nCol(); ++k)
            {
                m += A(i,k) * P(j,k);
            }
            chisq += GAPS_SQ((PSampler.mDMatrix(i,j) - m) / PSampler.mSMatrix(i,j));
        }
    }
    return chisq;
}

float GapsStatistics::meanChiSq(const SparseGibbsSampler &PSampler) const
{
    Matrix A(mAMeanMatrix / static_cast<float>(mStatUpdates));
    Matrix P(mPMeanMatrix / static_cast<float>(mStatUpdates));
    return 0.f; // TODO
}

Archive& operator<<(Archive &ar, GapsStatistics &stat)
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

