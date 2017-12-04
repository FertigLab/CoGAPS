#ifndef __COGAPS_MATRIX_ALGORITHM_H__
#define __COGAPS_MATRIX_ALGORITHM_H__

#include "Matrix.h"

namespace MatAlgo
{
    matrix_data_t sum(const Matrix &mat);
    matrix_data_t mean(const Matrix &mat);
    matrix_data_t nonZeroMean(const Matrix &mat);

    void matrixMultiplication(Matrix &result, const Matrix &A, const Matrix &B);
    void matrixSubtraction(Matrix &result, const Matrix &A, const Matrix &B);
    void elementWiseMatrixMultiplication(Matrix &result, const Matrix &A, const Matrix &B);
}

#endif