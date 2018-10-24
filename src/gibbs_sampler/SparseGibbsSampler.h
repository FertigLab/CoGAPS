#ifndef __COGAPS_SPARSE_GIBBS_SAMPLER_H__
#define __COGAPS_SPARSE_GIBBS_SAMPLER_H__

#include "GibbsSampler.h"

#include "../data_structures/HybridMatrix.h"
#include "../data_structures/SparseMatrix.h"
#include "../data_structures/SparseIterator.h"

#include <vector>

class GapsStatistics;

class SparseGibbsSampler : public GibbsSampler<SparseGibbsSampler, SparseMatrix, HybridMatrix>
{
public:

    friend class GapsStatistics;
    friend class GibbsSampler; // so impl()-> can access private members

    template <class DataType>
    SparseGibbsSampler(const DataType &data, bool transpose, bool subsetRows,
        float alpha, float maxGibbsMass, const GapsParameters &params,
        GapsRandomState *randState);

    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params);

    float chiSq() const;
    void sync(const SparseGibbsSampler &sampler, unsigned nThreads=1);
    void extraInitialization();

    friend Archive& operator<<(Archive &ar, const SparseGibbsSampler &s);
    friend Archive& operator>>(Archive &ar, SparseGibbsSampler &s);

#ifdef GAPS_DEBUG
    bool internallyConsistent() const;
#endif

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    Vector mZ1;
    Matrix mZ2;

    float mBeta;

    void changeMatrix(unsigned row, unsigned col, float delta);
    void safelyChangeMatrix(unsigned row, unsigned col, float delta);
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);

    void generateLookupTables();

    SparseGibbsSampler(const SparseGibbsSampler&); // = delete
    SparseGibbsSampler& operator=(const SparseGibbsSampler&); // = delete
};

template <class DataType>
SparseGibbsSampler::SparseGibbsSampler(const DataType &data, bool transpose,
bool subsetRows, float alpha, float maxGibbsMass, const GapsParameters &params,
GapsRandomState *randState)
    :
GibbsSampler(data, transpose, subsetRows, alpha, maxGibbsMass, params, randState),
mZ1(params.nPatterns),
mZ2(params.nPatterns, params.nPatterns),
mBeta(100.f)
{
    // check data for values less than 1
    for (unsigned j = 0; j < mDMatrix.nCol(); ++j)
    {
        SparseIterator it(mDMatrix.getCol(j));
        while (!it.atEnd())
        {
            if (it.getValue() < 1.f)
            {
                gaps_printf("\nError: Non-zero values less than 1 detected\n");
                gaps_printf("\n       Not allowed when useSparseOptimization is enabled\n");
                gaps_stop();
            }
            it.next();
        }
    }
}

// required function for the GibbsSampler interface
template <class DataType>
void SparseGibbsSampler::setUncertainty(const DataType &data, bool transpose,
bool subsetRows, const GapsParameters &params)
{
    // nop - SparseGibbsSampler assumes default uncertainty always
}

#endif // __COGAPS_SPARSE_GIIBS_SAMPLER_H__