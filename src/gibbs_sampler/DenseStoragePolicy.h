#ifndef __COGAPS_DENSE_STORAGE_POLICY_H__
#define __COGAPS_DENSE_STORAGE_POLICY_H__

#include "AlphaParameters.h"
#include "../GapsParameters.h"
#include "../data_structures/Matrix.h"
#include "../math/MatrixMath.h"
#include "../utils/Archive.h"

class DenseStorage
{
public:

    template <class DataType>
    DenseStorage(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params);

    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params);

    float chiSq() const;
    void sync(const DenseStorage &sampler, unsigned nThreads=1);
    void extraInitialization();
    
    float apSum() const;

    friend Archive& operator<<(Archive &ar, const DenseStorage &s);
    friend Archive& operator>>(Archive &ar, DenseStorage &s);

protected:

    void changeMatrix(unsigned row, unsigned col, float delta);
    void safelyChangeMatrix(unsigned row, unsigned col, float delta);
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    DenseStorage(const DenseStorage&); // = delete (no c++11)
    DenseStorage& operator=(const DenseStorage&); // = delete

    Matrix mDMatrix; // samples by genes for A, genes by samples for P
    Matrix mMatrix; // genes by patterns for A, samples by patterns for P
    const Matrix *mOtherMatrix; // pointer to P if this is A, and vice versa

    Matrix mSMatrix; // uncertainty values for each data point
    Matrix mAPMatrix; // cached product of A and P
};


template <class DataType>
DenseStorage::DenseStorage(const DataType &data, bool transpose,
bool subsetRows, const GapsParameters &params)
    :
mDMatrix(data, transpose, subsetRows, params.dataIndicesSubset),
mMatrix(mDMatrix.nCol(), params.nPatterns),
mOtherMatrix(NULL),
mSMatrix(gaps::pmax(mDMatrix, 0.1f)),
mAPMatrix(mDMatrix.nRow(), mDMatrix.nCol())
{
    mSMatrix.pad(1.f); // can't divide by zero
}

template <class DataType>
void DenseStorage::setUncertainty(const DataType &unc, bool transpose,
bool subsetRows, const GapsParameters &params)
{
    mSMatrix = Matrix(unc, transpose, subsetRows, params.dataIndicesSubset);
    mSMatrix.pad(1.f); // so that SIMD operations don't divide by zero
}

#endif // __COGAPS_DENSE_STORAGE_POLICY_H__