#ifndef __COGAPS_GAPS_STATISTICS_H__
#define __COGAPS_GAPS_STATISTICS_H__

#include "math/Math.h"
#include "math/MatrixMath.h"
#include "math/VectorMath.h"
#include "data_structures/Matrix.h"
#include "gibbs_sampler/DenseGibbsSampler.h"
#include "gibbs_sampler/SparseGibbsSampler.h"

#define GAPS_SQ(x) ((x) * (x))

class GapsStatistics
{
private:

    Matrix mAMeanMatrix;
    Matrix mAStdMatrix;
    Matrix mPMeanMatrix;
    Matrix mPStdMatrix;
    
    unsigned mStatUpdates;
    unsigned mNumPatterns;

public:

    GapsStatistics(unsigned nGenes, unsigned nSamples, unsigned nPatterns);

    template <class Sampler>
    void update(const Sampler &ASampler, const Sampler &PSampler);

    Matrix Amean() const;
    Matrix Pmean() const;
    Matrix Asd() const;
    Matrix Psd() const;
    
    float meanChiSq(const DenseGibbsSampler &PSampler) const;
    float meanChiSq(const SparseGibbsSampler &PSampler) const;

    // serialization
    friend Archive& operator<<(Archive &ar, const GapsStatistics &stat);
    friend Archive& operator>>(Archive &ar, GapsStatistics &stat);
};

template <class Sampler>
void GapsStatistics::update(const Sampler &ASampler, const Sampler &PSampler)
{
    mStatUpdates++;

    // update     
    // precision loss? use double?
    DEBUG_PING
    GAPS_ASSERT(mNumPatterns == ASampler.mMatrix.nCol());
    GAPS_ASSERT(mNumPatterns == PSampler.mMatrix.nCol());

    for (unsigned j = 0; j < mNumPatterns; ++j)
    {
        float norm = gaps::max(PSampler.mMatrix.getCol(j));
        norm = norm == 0.f ? 1.f : norm;
        GAPS_ASSERT(norm > 0.f);

        DEBUG_PING
        Vector quot(PSampler.mMatrix.getCol(j) / norm);
        GAPS_ASSERT(gaps::min(quot) >= 0.f);
        mPMeanMatrix.getCol(j) += quot;
        mPStdMatrix.getCol(j) += gaps::elementSq(quot);

        DEBUG_PING
        Vector prod(ASampler.mMatrix.getCol(j) * norm);
        GAPS_ASSERT(gaps::min(prod) >= 0.f);
        mAMeanMatrix.getCol(j) += prod;
        mAStdMatrix.getCol(j) += gaps::elementSq(prod);
    }
    DEBUG_PING
}

#endif