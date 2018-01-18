#ifndef __COGAPS_ALGORITHMS_H__
#define __COGAPS_ALGORITHMS_H__

#include "Matrix.h"

struct AlphaParameters
{
    float s;
    float su;
    
    AlphaParameters(float inS, float inSU)
        : s(inS), su(inSU)
    {}

    AlphaParameters& operator+(const AlphaParameters &other)
    {
        s += other.s;
        su -= other.su; // weird
        return *this;
    }
};

namespace gaps
{

namespace algo
{
    // vector algorithms    
    float sum(const Vector &vec);
    Vector scalarMultiple(const Vector &vec, float n);
    Vector squaredScalarMultiple(const Vector &vec, float n);
    Vector scalarDivision(const Vector &vec, float n);
    Vector squaredScalarDivision(const Vector &vec, float n);
    
    // generic matrix algorithms
    template<class GenericMatrix>
    float sum(const GenericMatrix &mat);

    template<class GenericMatrix>
    float mean(const GenericMatrix &mat);

    template<class GenericMatrix>
    float nonZeroMean(const GenericMatrix &mat);

    template<class GenericMatrix>
    GenericMatrix scalarDivision(const GenericMatrix &mat, float n);

    template<class GenericMatrix>
    GenericMatrix computeStdDev(const GenericMatrix &stdMat,
        const GenericMatrix &meanMat, unsigned nUpdates);

    // specific matrix algorithms
    bool isRowZero(const RowMatrix &mat, unsigned row);
    bool isColZero(const ColMatrix &mat, unsigned col);
    void matrixMultiplication(TwoWayMatrix &C, const ColMatrix &A,
        const RowMatrix &B);

    // chi2
    float loglikelihood(const TwoWayMatrix &D, const TwoWayMatrix &S,
        const TwoWayMatrix &AP);

    // change in likelihood
    float deltaLL(const MatrixChange &ch, const TwoWayMatrix &D,
        const TwoWayMatrix &S, const ColMatrix &A,
        const RowMatrix &P, const TwoWayMatrix &AP);

    // alpha parameters used in exchange and gibbsMass calculation
    AlphaParameters alphaParameters(const MatrixChange &ch,
        const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
        const RowMatrix &P, const TwoWayMatrix &AP);
}

}

template<class GenericMatrix>
float gaps::algo::sum(const GenericMatrix &mat)
{
    float sum = 0.f;
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            sum += mat(i,j);
        }
    }
    return sum;
}

template<class GenericMatrix>
float gaps::algo::mean(const GenericMatrix &mat)
{
    return gaps::algo::sum(mat) / (mat.nRow() * mat.nCol());
}

template<class GenericMatrix>
float gaps::algo::nonZeroMean(const GenericMatrix &mat)
{
    float sum = 0.f;
    unsigned nNonZeros = 0;
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            if (mat(i,j) != 0.f)
            {
                nNonZeros++;
                sum += mat(i,j);
            }
        }
    }
    return sum / (float)nNonZeros;
}

template<class GenericMatrix>
GenericMatrix gaps::algo::scalarDivision(const GenericMatrix &mat, float n)
{
    GenericMatrix retMat(mat.nRow(), mat.nCol());
    for (unsigned r = 0; r < retMat.nRow(); ++r)
    {
        for (unsigned c = 0; c < retMat.nCol(); ++c)
        {
            retMat(r,c) = mat(r,c) / n;
        }
    }
    return retMat;
}

template<class GenericMatrix>
GenericMatrix gaps::algo::computeStdDev(const GenericMatrix &stdMat,
const GenericMatrix &meanMat, unsigned nUpdates)
{
    GenericMatrix retMat(stdMat.nRow(), stdMat.nCol());
    float meanTerm = 0.0;
    for (unsigned r = 0; r < retMat.nRow(); ++r)
    {
        for (unsigned c = 0; c < retMat.nCol(); ++c)
        {
            meanTerm = meanMat(r,c) * meanMat(r,c) / (float)nUpdates;
            retMat(r,c) = std::sqrt((stdMat(r,c) - meanTerm) / ((float)nUpdates - 1.0));
        }
    }
    return retMat;
}

#endif