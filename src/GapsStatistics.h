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
public:

    GapsStatistics(unsigned nGenes, unsigned nSamples, unsigned nPatterns);

    template <class Sampler>
    void update(const Sampler &ASampler, const Sampler &PSampler);

    Matrix Amean() const;
    Matrix Pmean() const;
    Matrix Asd() const;
    Matrix Psd() const;

    void addChiSq(float chisq);
    void addAtomCount(unsigned atomA, unsigned atomP);

    std::vector<float> chisqHistory() const;
    std::vector<unsigned> atomHistory(char m) const;
    
    float meanChiSq(const DenseGibbsSampler &PSampler) const;
    float meanChiSq(const SparseGibbsSampler &PSampler) const;

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
    
    std::vector<float> mChisqHistory;
    std::vector<unsigned> mAtomHistoryA;
    std::vector<unsigned> mAtomHistoryP;

    unsigned mStatUpdates;
    unsigned mNumPatterns;
};

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