#include "VectorMath.h"
#include "MatrixMath.h"
#include "Math.h"
#include "../data_structures/SparseIterator.h"

float gaps::sparsity(const Matrix &mat)
{
    unsigned nNonZeroes = 0;
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            if (mat(i,j) > 0.f)
            {
                ++nNonZeroes;
            }
        }
    }
    float size = mat.nRow() * mat.nCol();
    return 1.f - static_cast<float>(nNonZeroes) / size;
}

float gaps::sparsity(const SparseMatrix &mat)
{
    unsigned nNonZeroes = 0;
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        SparseIterator<1> it(mat.getCol(j));
        while (!it.atEnd())
        {
            ++nNonZeroes;
            it.next();
        }
    }
    float size = mat.nRow() * mat.nCol();
    return 1.f - static_cast<float>(nNonZeroes) / size;
}

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

Matrix gaps::pmax(const Matrix & mat, float factor, float min_threshold)
{
    Matrix rmat(mat.nRow(), mat.nCol());
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            rmat(i,j) = gaps::max(mat(i,j) * factor, min_threshold);
        }   
    }
    return rmat;
}

//overload threshold=factor for back compatibility
Matrix gaps::pmax(const Matrix & mat, float factor)
{
   return gaps::pmax(mat, factor, factor);
}


Matrix operator*(const Matrix & mat, float f)
{
    Matrix rmat(mat.nRow(), mat.nCol());
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            rmat(i,j) = f * mat(i,j);
        }   
    }
    return rmat;
}

Matrix operator/(const Matrix & mat, float f)
{
    Matrix rmat(mat.nRow(), mat.nCol());
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            rmat(i,j) = mat(i,j) / f;
        }   
    }
    return rmat;
}
