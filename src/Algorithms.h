#ifndef __COGAPS_ALGORITHMS_H__
#define __COGAPS_ALGORITHMS_H__

#include "Matrix.h"

#if (!defined(GAPS_DEBUG) && defined(GAPS_INTERNAL_TESTS))
    #define GAPS_DEBUG
#endif

struct AlphaParameters
{
    double s;
    double su;
    
    AlphaParameters(double inS, double inSU)
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
    matrix_data_t sum(const Vector &vec);
    matrix_data_t sum(const TwoWayMatrix &mat);
    matrix_data_t sum(const RowMatrix &mat);
    matrix_data_t sum(const ColMatrix &mat);

    matrix_data_t mean(const TwoWayMatrix &mat);
    matrix_data_t nonZeroMean(const TwoWayMatrix &mat);

    ColMatrix computeStdDev(const ColMatrix &stdMat, const ColMatrix &meanMat,
        unsigned nUpdates);

    RowMatrix computeStdDev(const RowMatrix &stdMat, const RowMatrix &meanMat,
        unsigned nUpdates);

    Vector scalarMultiple(const Vector &vec, matrix_data_t n);
    Vector squaredScalarMultiple(const Vector &vec, matrix_data_t n);
    Vector scalarDivision(const Vector &vec, matrix_data_t n);
    Vector squaredScalarDivision(const Vector &vec, matrix_data_t n);

    ColMatrix scalarDivision(const ColMatrix &mat, matrix_data_t n);
    RowMatrix scalarDivision(const RowMatrix &mat, matrix_data_t n);

    bool isRowZero(const RowMatrix &mat, unsigned row);
    bool isColZero(const ColMatrix &mat, unsigned col);

    matrix_data_t loglikelihood(const TwoWayMatrix &D, const TwoWayMatrix &S,
        const TwoWayMatrix &AP);

    double deltaLL(const MatrixChange &ch, const TwoWayMatrix &D,
        const TwoWayMatrix &S, const ColMatrix &A,
        const RowMatrix &P, const TwoWayMatrix &AP);

    AlphaParameters alphaParameters(const MatrixChange &ch,
        const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
        const RowMatrix &P, const TwoWayMatrix &AP);

#ifdef GAPS_DEBUG
    bool checkAPMatrix(const TwoWayMatrix &AP, const ColMatrix &A,
        const RowMatrix &P);
#endif

}

}

#endif