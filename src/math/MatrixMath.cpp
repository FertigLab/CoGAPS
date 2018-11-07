#include "VectorMath.h"
#include "MatrixMath.h"
#include "Math.h"
#include "../data_structures/SparseIterator.h"

float gaps::nonZeroMean(const Matrix &mat)
{
    float sum = 0.f;
    unsigned nNonZeroes = 0;
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            sum += mat(i,j);
            if (mat(i,j) > 0.f)
            {
                ++nNonZeroes;
            }
        }
    }
    return sum / static_cast<float>(nNonZeroes);
}

float gaps::nonZeroMean(const SparseMatrix &mat)
{
    float sum = 0.f;
    unsigned nNonZeroes = 0;
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        SparseIterator<1> it(mat.getCol(j));
        while (!it.atEnd())
        {
            sum += get<1>(it);
            ++nNonZeroes;
            it.next();
        }
    }
    return sum / static_cast<float>(nNonZeroes);
}

Matrix gaps::pmax(Matrix mat, float p)
{
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            mat(i,j) = gaps::max(mat(i,j) * p, p);
        }   
    }
    return mat;
}

Matrix operator*(Matrix mat, float f)
{
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            mat(i,j) *= f;
        }   
    }
    return mat;
}

Matrix operator/(Matrix mat, float f)
{
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            mat(i,j) /= f;
        }   
    }
    return mat;
}