#include "Matrix.h"

#include <stdexcept>

/******************************** VECTOR *******************************/

double Vector::dotProduct(const Vector &vec) const
{
    if (size() != vec.size())
    {
        throw std::invalid_argument("vector dimensions not compatible");
    }

    double sum = 0.f;
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
        sum += mValues[i] * vec(i);
    }
    return sum;
}

double Vector::sum() const
{
    double sum = 0.f;
    for (unsigned i = 0; i < mValues.size(); ++i)
    {
        sum += mValues[i];
    }
    return sum;
}

/****************************** ROW MATRIX *****************************/

RowMatrix::RowMatrix(unsigned int nrow, unsigned int ncol)
: mNumRows(nrow), mNumCols(ncol)
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
    }
}

RowMatrix::RowMatrix(const Rcpp::NumericMatrix &rmat)
: mNumRows(rmat.nrow()), mNumCols(rmat.ncol())
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            mRows[i](j) = rmat(i,j);
        }
    }
}

RowMatrix::RowMatrix(const std::vector< std::vector<double> > &mat)
: mNumRows(mat.size()), mNumCols(mat[0].size())
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            mRows[i](j) = mat[i][j];
        }
    }
}

Vector& RowMatrix::getRow(unsigned int row)
{
    return mRows[row];
}

const Vector& RowMatrix::getRow(unsigned int row) const
{
    return mRows[row];
}

void RowMatrix::elemUpdate(const std::vector<ElementChange> &changes)
{
    for (unsigned i = 0; i < changes.size(); ++i)
    {
        mRows[changes[i].row](changes[i].col) += changes[i].delta;
    }
}

/**************************** COLUMN MATRIX ****************************/

ColMatrix::ColMatrix(unsigned int nrow, unsigned int ncol)
: mNumRows(nrow), mNumCols(ncol)
{
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        mCols.push_back(Vector(mNumRows));
    }
}

ColMatrix::ColMatrix(const Rcpp::NumericMatrix &rmat)
: mNumRows(rmat.nrow()), mNumCols(rmat.ncol())
{
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(Vector(mNumRows));
        for (unsigned i = 0; i < mNumRows; ++i)
        {
            mCols[i](j) = rmat(i,j);
        }
    }
}

Vector& ColMatrix::getCol(unsigned int col)
{
    return mCols[col];
}

void ColMatrix::elemUpdate(const std::vector<ElementChange> &changes)
{
    for (unsigned i = 0; i < changes.size(); ++i)
    {
        mCols[changes[i].col](changes[i].row) += changes[i].delta;
    }
}
