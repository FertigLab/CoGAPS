#include "Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/MatrixElement.h"
#include "../GapsAssert.h"

template <class MatA, class MatB>
inline void copyMatrix(MatA &dest, const MatB &source)
{
    GAPS_ASSERT(dest.nRow() == source.nRow());
    GAPS_ASSERT(dest.nCol() == source.nCol());

    for (unsigned i = 0; i < dest.nRow(); ++i)
    {
        for (unsigned j = 0; j < dest.nCol(); ++j)
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

RowMatrix::RowMatrix(const Matrix &mat, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
    : GenericMatrix(mat, transpose, partitionRows, indices)
{}

RowMatrix::RowMatrix(const std::string &path, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
    : GenericMatrix(path, transpose, partitionRows, indices)
{}

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

RowMatrix::RowMatrix(const ColMatrix &mat) : GenericMatrix(mat.nRow(), mat.nCol())
{
    copyMatrix(*this, mat);
}

RowMatrix& RowMatrix::operator=(const RowMatrix &mat)
{
    copyMatrix(*this, mat);
}

RowMatrix& RowMatrix::operator=(const ColMatrix &mat)
{
    copyMatrix(*this, mat);
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

ColMatrix::ColMatrix(const Matrix &mat, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
    : GenericMatrix(mat, transpose, partitionRows, indices)
{}

ColMatrix::ColMatrix(const std::string &path, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
    : GenericMatrix(path, transpose, partitionRows, indices)
{}

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

ColMatrix::ColMatrix(const RowMatrix &mat) : GenericMatrix(mat.nRow(), mat.nCol())
{
    copyMatrix(*this, mat);
}

ColMatrix& ColMatrix::operator=(const ColMatrix &mat)
{
    copyMatrix(*this, mat);
}

ColMatrix& ColMatrix::operator=(const RowMatrix &mat)
{
    copyMatrix(*this, mat);
}