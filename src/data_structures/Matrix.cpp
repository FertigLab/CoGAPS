#include "Matrix.h"

/******************************** HELPER *******************************/

template<class GenericMatrix>
static Rcpp::NumericMatrix convertToRMatrix(const GenericMatrix &mat)
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

template<class MatA, class MatB>
static void copyMatrix(MatA &dest, const MatB &source)
{
    for (unsigned i = 0; i < source.nRow(); ++i)
    {
        for (unsigned j = 0; j < source.nCol(); ++j)
        {
            dest(i,j) = source(i,j);
        }
    }
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
            this->operator()(i,j) = rmat(i,j);
        }
    }
}

void RowMatrix::operator=(const RowMatrix &mat)
{
    copyMatrix(*this, mat);
}

void RowMatrix::operator=(const ColMatrix &mat)
{
    copyMatrix(*this, mat);
}

Rcpp::NumericMatrix RowMatrix::rMatrix() const
{
    return convertToRMatrix(*this);
}

RowMatrix RowMatrix::operator/(float val) const
{
    RowMatrix mat(mNumRows, mNumCols);
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mat.getRow(i) = mRows[i] / val;
    }
    return mat;
}

Archive& operator<<(Archive &ar, RowMatrix &mat)
{
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        ar << mat.mRows[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, RowMatrix &mat)
{
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        ar >> mat.mRows[i];
    }
    return ar;
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
            this->operator()(i,j) = rmat(i,j);
        }
    }
}

ColMatrix ColMatrix::operator/(float val) const
{
    ColMatrix mat(mNumRows, mNumCols);
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mat.getCol(j) = mCols[j] / val;
    }
    return mat;
}

void ColMatrix::operator=(const ColMatrix &mat)
{
    copyMatrix(*this, mat);
}

void ColMatrix::operator=(const RowMatrix &mat)
{
    copyMatrix(*this, mat);
}

Rcpp::NumericMatrix ColMatrix::rMatrix() const
{
    return convertToRMatrix(*this);
}

Archive& operator<<(Archive &ar, ColMatrix &mat)
{
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        ar << mat.mCols[j];
    }
    return ar;
}

Archive& operator>>(Archive &ar, ColMatrix &mat)
{
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        ar >> mat.mCols[j];
    }
    return ar;
}