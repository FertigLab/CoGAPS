#include "Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/MatrixElement.h"
#include "../GapsAssert.h"

template <class MatA, class MatB>
void copyMatrix(MatA &dest, const MatB &source)
{
    GAPS_ASSERT(mNumRows == mat.nRow());
    GAPS_ASSERT(mNumCols == mat.nCol());

    for (unsigned i = 0; i < mNumRows; ++i)
    {
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            dest(i,j) = source(i,j);
        }
    }
}

/********************************** ROW MATRIX ********************************/

void RowMatrix::allocate()
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mData.push_back(Vector(mNumCols));
    }
}

RowMatrix::RowMatrix(unsigned nrow, unsigned ncol)
    : GenericMatrix(nrow, ncol)
{}

RowMatrix::RowMatrix(const Matrix &mat, bool transpose)
    : GenericMatrix(mat, transpose)
{}

RowMatrix::RowMatrix(const std::string &path, bool transpose)
    : GenericMatrix(path, transpose)
{}

RowMatrix::RowMatrix(const Matrix &mat, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
    : GenericMatrix(mat, transpose, partitionRows, indices)
{}

RowMatrix::RowMatrix(const std::string &path, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
    : GenericMatrix(path, transpose, partitionRows, indices)
{}

RowMatrix& RowMatrix::operator=(const ColMatrix &mat)
{
    copyMatrix(*this, mat);
}

float& RowMatrix::operator()(unsigned r, unsigned c)
{
    return mData[r][c];
}

float RowMatrix::operator()(unsigned r, unsigned c) const
{
    return mData[r][c];
}

Vector& RowMatrix::getRow(unsigned row)
{
    return mData[row];
}

const Vector& RowMatrix::getRow(unsigned row) const
{
    return mData[row];
}

float* RowMatrix::rowPtr(unsigned row)
{
    return mData[row].ptr();
}

const float* RowMatrix::rowPtr(unsigned row) const
{
    return mData[row].ptr();
}

/******************************** COLUMN MATRIX *******************************/

void ColMatrix::allocate()
{
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        mData.push_back(Vector(mNumRows));
    }
}

ColMatrix::ColMatrix(unsigned nrow, unsigned ncol)
    : GenericMatrix(nrow, ncol)
{}

ColMatrix::ColMatrix(const Matrix &mat, bool transpose)
    : GenericMatrix(mat, transpose)
{}

ColMatrix::ColMatrix(const std::string &path, bool transpose)
    : GenericMatrix(path, transpose)
{}

ColMatrix::ColMatrix(const Matrix &mat, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
    : GenericMatrix(mat, transpose, partitionRows, indices)
{}

ColMatrix::ColMatrix(const std::string &path, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
    : GenericMatrix(path, transpose, partitionRows, indices)
{}

ColMatrix& ColMatrix::operator=(const ColMatrix &mat)
{
    copyMatrix(*this, mat);
}

float& ColMatrix::operator()(unsigned r, unsigned c)
{
    return mData[c][r];
}

float ColMatrix::operator()(unsigned r, unsigned c) const
{
    return mData[c][r];
}

Vector& ColMatrix::getCol(unsigned col)
{
    return mData[col];
}

const Vector& ColMatrix::getCol(unsigned col) const
{
    return mData[col];
}

float* ColMatrix::colPtr(unsigned col)
{
    return mData[col].ptr();
}

const float* ColMatrix::colPtr(unsigned col) const
{
    return mData[col].ptr();
}
