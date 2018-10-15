#ifndef __COGAPS_MATRIX_MATH_H__
#define __COGAPS_MATRIX_MATH_H__

#include "../data_structures/Matrix.h"
#include "../data_structures/SparseMatrix.h"

namespace gaps
{
    float sum(const Matrix &mat);
    float sum(const SparseMatrix &mat);

    float mean(const Matrix &mat);
    float mean(const SparseMatrix &mat);
    
    float nonZeroMean(const Matrix &mat);
    float nonZeroMean(const SparseMatrix &mat);

    Matrix pmax(Matrix mat, float p);
}

Matrix operator*(Matrix mat, float f);
Matrix operator/(Matrix mat, float f);

#endif