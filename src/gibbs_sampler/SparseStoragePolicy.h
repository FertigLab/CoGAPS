#ifndef __COGAPS_SPARSE_STORAGE_POLICY_H__
#define __COGAPS_SPARSE_STORAGE_POLICY_H__

#include "AlphaParameters.h"
#include "../GapsParameters.h"
#include "../utils/Archive.h"
#include "../data_structures/HybridMatrix.h"
#include "../data_structures/SparseMatrix.h"

class SparseStorage
{
public:

    template <class DataType>
    SparseStorage(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params);

    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params);

    float chiSq() const;
    void sync(const SparseStorage &sampler, unsigned nThreads=1);
    void extraInitialization();

    float apSum() const;

    friend Archive& operator<<(Archive &ar, const SparseStorage &s);
    friend Archive& operator>>(Archive &ar, SparseStorage &s);

protected:

    void changeMatrix(unsigned row, unsigned col, float delta);
    void safelyChangeMatrix(unsigned row, unsigned col, float delta);
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);

    void generateLookupTables();

    SparseStorage(const SparseStorage&); // = delete
    SparseStorage& operator=(const SparseStorage&); // = delete

    SparseMatrix mDMatrix; // samples by genes for A, genes by samples for P
    HybridMatrix mMatrix; // genes by patterns for A, samples by patterns for P
    const HybridMatrix *mOtherMatrix; // pointer to P if this is A, and vice versa

    Matrix mZ2;
    Vector mZ1;

    float mBeta;
};

template <class DataType>
SparseStorage::SparseStorage(const DataType &data, bool transpose,
bool subsetRows, const GapsParameters &params)
    :
mDMatrix(data, transpose, subsetRows, params.dataIndicesSubset),
mMatrix(mDMatrix.nCol(), params.nPatterns),
mOtherMatrix(NULL),
mZ2(params.nPatterns, params.nPatterns),
mZ1(params.nPatterns),
mBeta(100.f)
{}

// required function for the GibbsSampler interface
template <class DataType>
void SparseStorage::setUncertainty(const DataType &data, bool transpose,
bool subsetRows, const GapsParameters &params)
{
    // nop - SparseGibbsSampler assumes default uncertainty always
}

#endif // __COGAPS_SPARSE_STORAGE_POLICY_H__