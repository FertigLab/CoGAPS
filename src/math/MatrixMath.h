#ifndef __COGAPS_MATRIX_MATH_H__
#define __COGAPS_MATRIX_MATH_H__

// in order for overload resolution to work correctly with vectors we need
// to include VectorMath, code won't compile if overload resolution fails
#include "VectorMath.h"

#include "../data_structures/Matrix.h"
#include "../data_structures/HybridMatrix.h"
#include "../data_structures/SparseMatrix.h"

namespace gaps
{
    template <class MatrixType>
    float min(const MatrixType &mat);

    template <class MatrixType>
    float max(const MatrixType &mat);

    template <class MatrixType>
    float sum(const MatrixType &mat);

    template <class MatrixType>
    float mean(const MatrixType &mat);
    
    float nonZeroMean(const Matrix &mat);
    float nonZeroMean(const SparseMatrix &mat);

    Matrix pmax(Matrix mat, float p);
} // namespace gaps

Matrix operator*(Matrix mat, float f);
Matrix operator/(Matrix mat, float f);

template <class MatrixType>
float gaps::min(const MatrixType &mat)
{
    float mn = 0.f;
    for (unsigned i = 0; i < mat.nCol(); ++i)
    {
        float cmin = gaps::min(mat.getCol(i));
        mn = (cmin < mn) ? cmin : mn;
    }
    return mn;
}

template <class MatrixType>
float gaps::max(const MatrixType &mat)
{
    float mx = 0.f;
    for (unsigned i = 0; i < mat.nCol(); ++i)
    {
        float cmax = gaps::max(mat.getCol(i));
        mx = (cmax < mx) ? cmax : mx;
    }
    return mx;
}

template <class MatrixType>
float gaps::sum(const MatrixType &mat)
{
    float sum = 0.f;
    for (unsigned i = 0; i < mat.nCol(); ++i)
    {
        sum += gaps::sum(mat.getCol(i));
    }
    return sum;
}

template <class MatrixType>
float gaps::mean(const MatrixType &mat)
{
    return gaps::sum(mat) / (mat.nRow() * mat.nCol());
}

#endif