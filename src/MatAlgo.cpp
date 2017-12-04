#include "MatAlgo.h"
#include "Matrix.h"

#include <algorithm>

matrix_data_t MatAlgo::mean(const RowMatrix &mat)
{
    return MatAlgo::sum(mat) / (mat.nRow() * mat.nCol());
}

matrix_data_t MatAlgo::sum(const RowMatrix &mat)
{
    matrix_data_t sum = 0.0;
    for (unsigned r = 0; r < mat.nRow(); ++r)
    {
        sum += mat.getRow(r).sum();
    }
    return sum;
}

matrix_data_t MatAlgo::nonZeroMean(const Matrix &mat)
{
    matrix_data_t sum = 0.0;
    unsigned int nNonZeros = 0;
    for (unsigned int r = 0; r < mat.nRow(); ++r)
    {
        for (unsigned int c = 0; c < mat.nCol(); ++c)
        {
            if (mat(r,c) != 0.0)
            {
                nNonZeros++;
                sum += mat(r,c);
            }
        }
    }
    return sum / nNonZeros;
}


void MatAlgo::matrixMultiplication(Matrix &result, const Matrix &A, const Matrix &B)
{
    if (result.nRow() != A.nRow() ||
        result.nCol() != B.nCol() ||
        A.nCol() != B.nRow())
    {
        throw std::invalid_argument("invalid matrix dimensions");
    }

    for (unsigned int r = 0; r < result.nRow(); ++r)
    {
        for (unsigned int c = 0; c < result.nCol(); ++c)
        {
            result(r,c) = 0;
            for (unsigned int k = 0; k < A.nCol(); ++k)
            {
                result(r,c) += A(r,k) * B(k,c);
            }
        }
    }
}

void MatAlgo::matrixSubtraction(Matrix &result, const Matrix &A, const Matrix &B)
{
    if (result.nRow() != A.nRow() ||
        result.nRow() != B.nRow() ||
        result.nCol() != A.nCol() ||
        result.nCol() != B.nCol() ||
        A.nRow() != B.nRow() ||
        A.nCol() != B.nCol())
    {
        throw std::invalid_argument("invalid matrix dimensions");
    }

    for (unsigned int r = 0; r < result.nRow(); ++r)
    {
        for (unsigned int c = 0; c < result.nCol(); ++c)
        {
            result(r,c) = A(r,c) - B(r,c);
        }
    }
}

void MatAlgo::elementWiseMatrixMultiplication(Matrix &result, const Matrix &A, const Matrix &B)
{
    if (result.nRow() != A.nRow() ||
        result.nRow() != B.nRow() ||
        result.nCol() != A.nCol() ||
        result.nCol() != B.nCol() ||
        A.nRow() != B.nRow() ||
        A.nCol() != B.nCol())
    {
        throw std::invalid_argument("invalid matrix dimensions");
    }

    for (unsigned int r = 0; r < result.nRow(); ++r)
    {
        for (unsigned int c = 0; c < result.nCol(); ++c)
        {
            result(r,c) = A(r,c) * B(r,c);
        }
    }
}
