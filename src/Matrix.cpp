#include "Matrix.h"

#include <stdexcept>

static const float EPSILON = 1.e-10;

/******************************** HELPER *******************************/

static void updateHelper2(float& val, float delta)
{
    val = std::abs(val + delta) < EPSILON ? 0.0 : val + delta;
}

template<class GenericMatrix>
static void updateHelper(GenericMatrix &mat, const MatrixChange &change)
{
    updateHelper2(mat(change.row1, change.col1), change.delta1);
    if (change.nChanges > 1)
    {
        updateHelper2(mat(change.row2, change.col2), change.delta2);
    }
}

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

/******************************** VECTOR *******************************/

void Vector::concat(const Vector& vec)
{
    mValues.insert(mValues.end(), vec.mValues.begin(), vec.mValues.end());
}

void Vector::operator+=(const Vector &vec)
{
    for (unsigned i = 0; i < size(); ++i)
    {
        mValues[i] += vec[i];
    }
}

void Vector::operator=(const Vector &vec)
{
    mValues = vec.mValues;
}

Vector Vector::operator+(Vector vec) const
{
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] += mValues[i];
    }
    return vec;
}

Vector Vector::operator-(Vector vec) const
{
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] -= mValues[i];
    }
    return vec;
}

Vector Vector::operator*(Vector vec) const
{
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] *= mValues[i];
    }
    return vec;
}

Vector Vector::operator/(Vector vec) const
{
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] /= mValues[i];
    }
    return vec;
}

Vector Vector::operator+(float val) const
{
    Vector vec(*this);
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] += val;
    }
    return vec;
}

Vector Vector::operator-(float val) const
{
    Vector vec(*this);
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] -= val;
    }
    return vec;
}

Vector Vector::operator*(float val) const
{
    Vector vec(*this);
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] *= val;
    }
    return vec;
}

Vector Vector::operator/(float val) const
{
    Vector vec(*this);
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] /= val;
    }
    return vec;
}

Archive& operator<<(Archive &ar, Vector &vec)
{
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        ar << vec[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, Vector &vec)
{
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        ar >> vec.mValues[i];
    }
    return ar;
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

RowMatrix::RowMatrix(const RowMatrix &mat)
: mNumRows(mat.nRow()), mNumCols(mat.nCol())
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(mat.getRow(i));
    }
}

RowMatrix::RowMatrix(const ColMatrix &mat)
: mNumRows(mat.nRow()), mNumCols(mat.nCol())
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            this->operator()(i,j) = mat(i,j);
        }
    }
}

void RowMatrix::update(const MatrixChange &change)
{
    updateHelper<RowMatrix>(*this, change);
}

Rcpp::NumericMatrix RowMatrix::rMatrix() const
{
    return convertToRMatrix<RowMatrix>(*this);
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

ColMatrix::ColMatrix(const ColMatrix &mat)
: mNumRows(mat.nRow()), mNumCols(mat.nCol())
{
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(mat.getCol(j));
    }
}

void ColMatrix::update(const MatrixChange &change)
{
    updateHelper<ColMatrix>(*this, change);
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

Rcpp::NumericMatrix ColMatrix::rMatrix() const
{
    return convertToRMatrix<ColMatrix>(*this);
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