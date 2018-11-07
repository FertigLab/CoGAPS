#include "VectorMath.h"
#include "MatrixMath.h"
#include "Math.h"
#include "../data_structures/SparseIterator.h"

float gaps::min(const Matrix &mat)
{
    float mn = mat(0,0);
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        mn = gaps::min(gaps::min(mat.getCol(j)), mn);
    }
    return mn;
}

float gaps::min(const SparseMatrix &mat)
{
    float mn = 0.f;
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        SparseIterator<1> it(mat.getCol(j));
        while (!it.atEnd())
        {
            mn = (get<1>(it) < mn) ? get<1>(it) : mn;
        }
    }
    return mn;
}

float gaps::max(const Matrix &mat)
{
    float mx = mat(0,0);
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        mx = gaps::max(gaps::max(mat.getCol(j)), mx);
    }
    return mx;
}

float gaps::max(const SparseMatrix &mat)
{
    float mx = 0.f;
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        SparseIterator<1> it(mat.getCol(j));
        while (!it.atEnd())
        {
            mx = (get<1>(it) > mx) ? get<1>(it) : mx;
            it.next();
        }
    }
    return mx;
}

float gaps::sum(const Matrix &mat)
{
    float sum = 0.f;
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        sum += gaps::sum(mat.getCol(j));
    }
    return sum;
}

float gaps::sum(const HybridMatrix &mat)
{
    float sum = 0.f;
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        sum += gaps::sum(mat.getCol(j));
    }
    return sum;
}

float gaps::sum(const SparseMatrix &mat)
{
    float sum = 0.f;
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        SparseIterator<1> it(mat.getCol(j));
        while (!it.atEnd())
        {
            sum += get<1>(it);
            it.next();
        }
    }
    return sum;
}

float gaps::mean(const Matrix &mat)
{
    return gaps::sum(mat) / (mat.nRow() * mat.nCol());
}

float gaps::mean(const SparseMatrix &mat)
{
    return gaps::sum(mat) / (mat.nRow() * mat.nCol());
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