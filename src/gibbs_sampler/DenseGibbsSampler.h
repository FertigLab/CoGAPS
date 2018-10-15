#ifndef __COGAPS_DENSE_GIBBS_SAMPLER_H__
#define __COGAPS_DENSE_GIBBS_SAMPLER_H__

#include "../data_structures/Matrix.h"
#include "GibbsSampler.h"

class GapsStatistics;

class DenseGibbsSampler : public GibbsSampler<DenseGibbsSampler, Matrix, Matrix>
{
public:

    friend class GapsStatistics;
    friend class GibbsSampler; // so impl()-> can access private members

    template <class DataType>
    DenseGibbsSampler(const DataType &data, bool transpose, bool subsetRows,
        float alpha, float maxGibbsMass, const GapsParameters &params);

    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, bool subsetRows,
        const GapsParameters &params);

    float chiSq() const;
    void sync(const DenseGibbsSampler &sampler, unsigned nThreads=1);
    void recalculateAPMatrix();

    friend Archive& operator<<(Archive &ar, DenseGibbsSampler &s);
    friend Archive& operator>>(Archive &ar, DenseGibbsSampler &s);

private:

    Matrix mSMatrix; // uncertainty values for each data point
    Matrix mAPMatrix; // cached product of A and P

    void changeMatrix(unsigned row, unsigned col, float delta);
    void safelyChangeMatrix(unsigned row, unsigned col, float delta);
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);

    void updateAPMatrix(unsigned row, unsigned col, float delta);

    DenseGibbsSampler(const DenseGibbsSampler&); // = delete (no c++11)
    DenseGibbsSampler& operator=(const DenseGibbsSampler&); // = delete
};

template <class DataType>
DenseGibbsSampler::DenseGibbsSampler(const DataType &data, bool transpose,
bool subsetRows, float alpha, float maxGibbsMass, const GapsParameters &params)
    :
GibbsSampler(data, transpose, subsetRows, alpha, maxGibbsMass, params),
mSMatrix(gaps::pmax(mDMatrix, 0.1f)),
mAPMatrix(mDMatrix.nRow(), mDMatrix.nCol())
{}

template <class DataType>
void DenseGibbsSampler::setUncertainty(const DataType &unc, bool transpose,
bool subsetRows, const GapsParameters &params)
{
    mSMatrix = Matrix(unc, transpose, subsetRows, params.dataIndicesSubset);
}

#endif // __COGAPS_DENSE_GIIBS_SAMPLER_H__