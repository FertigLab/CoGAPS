#include "Matrix.h"

#include <stdexcept>

static const double EPSILON = 1.e-10;

/******************************** HELPER *******************************/

static void updateHelper(double& val, double delta)
{
    val += delta;
    if (std::abs(val) < EPSILON)
    {
        val = 0.0;
    }
}

template <class GenericMatrix>
Rcpp::NumericMatrix convertToRMatrix(const GenericMatrix &mat)
{
    Rcpp::NumericMatrix rmat(mat.nRow(), mat.nCol());
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            rmat(i,j) = mat(i,j);
        }
    }
    return rmat;
}

/******************************** VECTOR *******************************/

Rcpp::NumericVector Vector::rVec() const
{
    return Rcpp::wrap(mValues);
}

void Vector::concat(const Vector& vec)
{
    mValues.insert(mValues.end(), vec.mValues.begin(), vec.mValues.end());
}

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

RowMatrix::RowMatrix(const RowMatrix &mat)
: mNumRows(mat.nRow()), mNumCols(mat.nCol())
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(mat.getRow(i));
    }
}

void RowMatrix::update(const MatrixChange &change)
{
    updateHelper(mRows[change.row1][change.col1], change.delta1);
    if (change.nChanges > 1)
    {
        updateHelper(mRows[change.row2][change.col2], change.delta2);
    }
}

Rcpp::NumericMatrix RowMatrix::rMatrix() const
{
    return convertToRMatrix<RowMatrix>(*this);
}

void Vector::operator+=(const Vector &vec)
{
    for (unsigned i = 0; i < size(); ++i)
    {
        mValues[i] += vec(i);
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

ColMatrix::ColMatrix(const ColMatrix &mat)
: mNumRows(mat.nRow()), mNumCols(mat.nCol())
{
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(mat.getCol(j));
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
    updateHelper(mCols[change.col1][change.row1], change.delta1);
    if (change.nChanges > 1)
    {
        updateHelper(mCols[change.col2][change.row2], change.delta2);
    }
}

Rcpp::NumericMatrix ColMatrix::rMatrix() const
{
    return convertToRMatrix<ColMatrix>(*this);
}
