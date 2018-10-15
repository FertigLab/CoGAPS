#ifndef __COGAPS_SPARSE_GIBBS_SAMPLER_H__
#define __COGAPS_SPARSE_GIBBS_SAMPLER_H__

#include "GibbsSampler.h"

#include "../data_structures/HybridMatrix.h"
#include "../data_structures/SparseMatrix.h"

#include <vector>

class GapsStatistics;

class SparseGibbsSampler : public GibbsSampler<SparseGibbsSampler, SparseMatrix, HybridMatrix>
{
public:

    friend class GapsStatistics;
    friend class GibbsSampler; // so impl()-> can access private members

    template <class DataType>
    SparseGibbsSampler(const DataType &data, bool transpose, bool subsetRows,
        float alpha, float maxGibbsMass, const GapsParameters &params);

    float chiSq() const;
    void sync(const SparseGibbsSampler &sampler, unsigned nThreads=1);

    friend Archive& operator<<(Archive &ar, SparseGibbsSampler &s);
    friend Archive& operator>>(Archive &ar, SparseGibbsSampler &s);

private :

    std::vector<float> mZ1;
    std::vector<float> mZ2;

    float mBeta;

    void changeMatrix(unsigned row, unsigned col, float delta);
    void safelyChangeMatrix(unsigned row, unsigned col, float delta);
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);

    float& Z1(unsigned pattern);
    float& Z2(unsigned pattern1, unsigned pattern2);
    void generateLookupTables();

    SparseGibbsSampler(const SparseGibbsSampler&); // = delete
    SparseGibbsSampler& operator=(const SparseGibbsSampler&); // = delete
};

template <class DataType>
SparseGibbsSampler::SparseGibbsSampler(const DataType &data, bool transpose,
bool subsetRows, float alpha, float maxGibbsMass, const GapsParameters &params)
    :
GibbsSampler(data, transpose, subsetRows, alpha, maxGibbsMass, params),
mBeta(100.f)
{}

#endif // __COGAPS_SPARSE_GIIBS_SAMPLER_H__