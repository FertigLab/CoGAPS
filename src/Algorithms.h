#ifndef __COGAPS_ALGORITHMS_H__
#define __COGAPS_ALGORITHMS_H__

#include "Matrix.h"
#include <xmmintrin.h>

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

class vec4f
{
private:
    
    __m128 mValue;

public:

    inline vec4f() {}
    inline vec4f(float f) : mValue(_mm_set1_ps(f)) {}
    inline vec4f(float f0, float f1, float f2, float f3)
        : mValue(_mm_setr_ps(f0,f1,f2,f3))
    {}
    inline vec4f(const __m128 &rhs) : mValue(rhs) {}
    
    inline vec4f& operator=(const __m128 &rhs)
    {
        mValue = rhs;
        return *this;
    }

    inline operator __m128() const {return mValue;}

    inline vec4f& operator+=(const vec4f &rhs)
    {
        *this = *this + rhs;
        return *this;
    }
};

inline vec4f operator+(const vec4f &lhs, const vec4f &rhs)
{
    return _mm_add_ps(lhs,rhs);
}

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

    void matrixMultiplication(TwoWayMatrix &C, const ColMatrix &A,
        const RowMatrix &B);

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

    float deltaLL(const MatrixChange &ch, const TwoWayMatrix &D,
        const TwoWayMatrix &S, const ColMatrix &A,
        const RowMatrix &P, const TwoWayMatrix &AP);

    AlphaParameters alphaParameters(const MatrixChange &ch,
        const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
        const RowMatrix &P, const TwoWayMatrix &AP);
}

}

#endif