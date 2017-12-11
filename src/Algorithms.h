#ifndef __COGAPS_ALGORITHMS_H__
#define __COGAPS_ALGORITHMS_H__

#include "Matrix.h"

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
    matrix_data_t mean(const TwoWayMatrix &mat);
    matrix_data_t nonZeroMean(const TwoWayMatrix &mat);
    bool isRowZero(const RowMatrix &mat, unsigned row);
    bool isColZero(const ColMatrix &mat, unsigned col);

    double deltaLL(const MatrixChange &ch, const TwoWayMatrix &D,
        const TwoWayMatrix &S, const ColMatrix &A,
        const RowMatrix &P, const TwoWayMatrix &AP);

    AlphaParameters alphaParameters(const MatrixChange &ch,
        const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
        const RowMatrix &P, const TwoWayMatrix &AP);
}

}

#endif