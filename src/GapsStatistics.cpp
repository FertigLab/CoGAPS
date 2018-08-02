#include "GapsStatistics.h"
#include "math/Algorithms.h"

GapsStatistics::GapsStatistics(unsigned nRow, unsigned nCol, unsigned nPatterns,
PumpThreshold t)
    :
mAMeanMatrix(nRow, nPatterns), mAStdMatrix(nRow, nPatterns),
mPMeanMatrix(nPatterns, nCol), mPStdMatrix(nPatterns, nCol),
mStatUpdates(0), mNumPatterns(nPatterns), mPumpMatrix(nRow, nCol),
mPumpThreshold(t), mPumpStatUpdates(0)
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

ColMatrix GapsStatistics::Amean() const
{
    return mAMeanMatrix / mStatUpdates;
}

ColMatrix GapsStatistics::Asd() const
{
    return gaps::algo::computeStdDev(mAStdMatrix, mAMeanMatrix,
        mStatUpdates);
}

RowMatrix GapsStatistics::Pmean() const
{
    return mPMeanMatrix / mStatUpdates;
}

RowMatrix GapsStatistics::Psd() const
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

static unsigned geneThreshold(const ColMatrix &rankMatrix, unsigned pat)
{
    float cutRank = rankMatrix.nRow();
    for (unsigned i = 0; i < rankMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < rankMatrix.nCol(); ++j)
        {
            if (j != pat && rankMatrix(i,j) <= rankMatrix(i,pat))
            {
                cutRank = gaps::min(cutRank, gaps::max(0.f, rankMatrix(i,pat)-1));
            }
        }
    }
    return static_cast<unsigned>(cutRank);
}

void GapsStatistics::patternMarkers(ColMatrix normedA, RowMatrix normedP,
ColMatrix &statMatrix)
{
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
            statMatrix(i,minNdx)++;
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
                    statMatrix(i,j)++;
                }
            }
        }
    }
}

void GapsStatistics::updatePump(const AmplitudeGibbsSampler &ASampler,
const PatternGibbsSampler &PSampler)
{
    mPumpStatUpdates++;
    patternMarkers(ASampler.mMatrix, PSampler.mMatrix, mPumpMatrix);
}

RowMatrix GapsStatistics::pumpMatrix() const
{
    unsigned denom = mPumpStatUpdates != 0 ? mPumpStatUpdates : 1.f;
    return RowMatrix(mPumpMatrix / denom);
}

RowMatrix GapsStatistics::meanPattern()
{
    ColMatrix Amean(mAMeanMatrix / static_cast<float>(mStatUpdates));
    RowMatrix Pmean(mPMeanMatrix / static_cast<float>(mStatUpdates));
    ColMatrix mat(mAMeanMatrix.nRow(), mAMeanMatrix.nCol());
    patternMarkers(Amean, Pmean, mat);
    return RowMatrix(mat);
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

