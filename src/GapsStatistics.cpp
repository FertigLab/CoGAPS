#include "GapsStatistics.h"
#include "math/Algorithms.h"

GapsStatistics::GapsStatistics(unsigned nRow, unsigned nCol, unsigned nFactor)
    : mAMeanMatrix(nRow, nFactor), mAStdMatrix(nRow, nFactor),
        mPMeanMatrix(nFactor, nCol), mPStdMatrix(nFactor, nCol),
        mStatUpdates(0), mNumPatterns(nFactor), mPumpMatrix(nRow, nCol),
        mPumpThreshold(PUMP_CUT), mPumpStatUpdates(0)
{}

void GapsStatistics::update(const AmplitudeGibbsSampler &ASampler,
const PatternGibbsSampler &PSampler)
{
    mStatUpdates++;

    // update     
    for (unsigned j = 0; j < mNumPatterns; ++j)
    {
        float norm = gaps::algo::sum(PSampler.mMatrix.getRow(j));
        norm = norm == 0.f ? 1.f : norm;

        Vector quot(PSampler.mMatrix.getRow(j) / norm);
        mPMeanMatrix.getRow(j) += quot;
        mPStdMatrix.getRow(j) += gaps::algo::elementSq(quot);

        Vector prod(ASampler.mMatrix.getCol(j) * norm);
        mAMeanMatrix.getCol(j) += prod;
        mAStdMatrix.getCol(j) += gaps::algo::elementSq(prod); 
    }
}

static unsigned geneThreshold(const ColMatrix &rankMatrix, unsigned pat)
{
    float cutRank = rankMatrix.nRow();
    for (unsigned i = 0; i < rankMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < rankMatrix.nCol(); ++j)
        {
            if (j != pat && rankMatrix(i,j) <= rankMatrix(i,pat))
            {
                cutRank = std::min(cutRank, std::max(0.f, rankMatrix(i,pat)-1));
            }
        }
    }
    return static_cast<unsigned>(cutRank);
}

void GapsStatistics::updatePump(const AmplitudeGibbsSampler &ASampler,
const PatternGibbsSampler &PSampler)
{
    mPumpStatUpdates++;

    // copy matrices so we can modify locally
    ColMatrix normedA(ASampler.mMatrix);
    RowMatrix normedP(PSampler.mMatrix);

    // helpful notation
    unsigned nGenes = normedA.nRow();
    unsigned nPatterns = normedA.nCol();

    // norm matrices
    for (unsigned j = 0; j < nPatterns; ++j)
    {
        float factor = gaps::algo::sum(normedP.getRow(j));
        factor = (factor == 0) ? 1.f : factor;

        normedP.getRow(j) = normedP.getRow(j) / factor;
        float scale = gaps::algo::max(normedP.getRow(j));
        normedA.getCol(j) = normedA.getCol(j) * factor * scale;
    }

    RowMatrix scaledA(nGenes, nPatterns);
    for (unsigned i = 0; i < nGenes; ++i)
    {
        for (unsigned j = 0; j < nPatterns; ++j)
        {
            scaledA(i,j) = normedA(i,j);
        }
    }
   
    // compute sstat
    RowMatrix rsStat(nGenes, nPatterns);
    ColMatrix csStat(nGenes, nPatterns);
    Vector lp(nPatterns), diff(nPatterns);
    for (unsigned j = 0; j < nPatterns; ++j)
    {
        lp[j] = 1.f;
        for (unsigned i = 0; i < nGenes; ++i)
        {
            float geneMax = gaps::algo::max(scaledA.getRow(i));
            diff = geneMax > 0.f ? scaledA.getRow(i) / geneMax - lp : lp * -1.f;
            rsStat(i,j) = std::sqrt(gaps::algo::dot(diff, diff));
            csStat(i,j) = rsStat(i,j);
        }
        lp[j] = 0.f;
    }

    // update PUMP matrix
    if (mPumpThreshold == PUMP_UNIQUE)
    {
        for (unsigned i = 0; i < nGenes; ++i)
        {
            unsigned minNdx = gaps::algo::whichMin(rsStat.getRow(i));
            mPumpMatrix(i,minNdx)++;
        }
    }
    else if (mPumpThreshold == PUMP_CUT)
    {
        ColMatrix rankMatrix(nGenes, nPatterns);
        for (unsigned j = 0; j < nPatterns; ++j)
        {
            rankMatrix.getCol(j) = gaps::algo::rank(csStat.getCol(j));
        }
        
        for (unsigned j = 0; j < nPatterns; ++j)
        {
            unsigned cutRank = geneThreshold(rankMatrix, j);
            for (unsigned i = 0; i < nGenes; ++i)
            {
                if (rankMatrix(i,j) <= cutRank)
                {
                    mPumpMatrix(i,j)++;
                }
            }
        }
    }
}

ColMatrix GapsStatistics::AMean() const
{
    return mAMeanMatrix / mStatUpdates;
}

ColMatrix GapsStatistics::AStd() const
{
    return gaps::algo::computeStdDev(mAStdMatrix, mAMeanMatrix,
        mStatUpdates);
}

RowMatrix GapsStatistics::PMean() const
{
    return mPMeanMatrix / mStatUpdates;
}

RowMatrix GapsStatistics::PStd() const
{
    return gaps::algo::computeStdDev(mPStdMatrix, mPMeanMatrix,
        mStatUpdates);
}

float GapsStatistics::meanChiSq(const AmplitudeGibbsSampler &ASampler) const
{
    ColMatrix A = mAMeanMatrix / mStatUpdates;
    RowMatrix P = mPMeanMatrix / mStatUpdates;
    RowMatrix M(gaps::algo::matrixMultiplication(A, P));
    return 2.f * gaps::algo::loglikelihood(ASampler.mDMatrix, ASampler.mSMatrix,
        M);
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

