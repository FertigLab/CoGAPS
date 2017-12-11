#include "Matrix.h"

#include <stdexcept>

/****************************** ROW MATRIX *****************************/

RowMatrix::RowMatrix(unsigned nrow, unsigned ncol)
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

Vector& RowMatrix::getRow(unsigned row)
{
    return mRows[row];
}

const Vector& RowMatrix::getRow(unsigned row) const
{
    return mRows[row];
}

void RowMatrix::update(const MatrixChange &change)
{
    mRows[change.row1](change.col1) += change.delta1;
    if (change.nChanges > 1)
    {
        mRows[change.row2](change.col2) += change.delta2;
    }
}

/**************************** COLUMN MATRIX ****************************/

ColMatrix::ColMatrix(unsigned nrow, unsigned ncol)
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
            mCols[j](i) = rmat(i,j);
        }
    }
}

Vector& ColMatrix::getCol(unsigned col)
{
    return mCols[col];
}

const Vector& ColMatrix::getCol(unsigned col) const
{
    return mCols[col];
}

void ColMatrix::update(const MatrixChange &change)
{
    mCols[change.col1](change.row1) += change.delta1;
    if (change.nChanges > 1)
    {
        mCols[change.col2](change.row2) += change.delta2;
    }
}
