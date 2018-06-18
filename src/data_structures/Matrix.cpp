#include "Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/MatrixElement.h"

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

// if partitionRows is false, partition columns instead
// rows of matrix should be partition dimension, i.e. need to transpose
// is partitionRows is false
template <class Matrix>
static void fill(Matrix &mat, FileParser &p, bool partitionRows, std::vector<unsigned> whichIndices)
{
    unsigned rowSelect = partitionRows ? 0 : 1;
    unsigned colSelect = partitionRows ? 1 : 0;

    while (p.hasNext())
    {
        MatrixElement e(p.getNext());
        std::vector<unsigned>::iterator newRowIndex = std::find(whichIndices.begin(), whichIndices.end(), e[rowSelect]);
        if (newRowIndex != whichIndices.end())
        {
            unsigned row = std::distance(whichIndices.begin(), newRowIndex);
            mat.operator()(row, e[colSelect]) = e.value();
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

RowMatrix::RowMatrix(FileParser &p) : mNumRows(p.nRow()), mNumCols(p.nCol())
{
    // allocate matrix
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
    }

    // populate matrix
    while (p.hasNext())
    {
        MatrixElement e(p.getNext());
        this->operator()(e.row(), e.col()) = e.value();
    }
}

RowMatrix::RowMatrix(FileParser &p, bool partitionRows, std::vector<unsigned> whichIndices)
{
    mNumRows = whichIndices.size();
    mNumCols = partitionRows ? p.nCol() : p.nRow();

    // allocate matrix
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
    }

    // fill in matrix
    fill(*this, p, partitionRows, whichIndices);
}

RowMatrix& RowMatrix::operator=(const ColMatrix &mat)
{
    copyMatrix(*this, mat);
    return *this;
}

RowMatrix RowMatrix::operator*(float val) const
{
    RowMatrix mat(mNumRows, mNumCols);
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mat.getRow(i) = mRows[i] * val;
    }
    return mat;
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

RowMatrix RowMatrix::pmax(float scale) const
{
    RowMatrix mat(mNumRows, mNumCols);
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            mat(i,j) = std::max(this->operator()(i,j) * scale, scale);
        }
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

ColMatrix::ColMatrix(const RowMatrix &mat)
: mNumRows(mat.nRow()), mNumCols(mat.nCol())
{
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(Vector(mNumRows));
        for (unsigned i = 0; i < mNumRows; ++i)
        {
            this->operator()(i,j) = mat(i,j);
        }
    }
}

ColMatrix::ColMatrix(FileParser &p) : mNumRows(p.nRow()), mNumCols(p.nCol())
{
    // allocate matrix
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(Vector(mNumRows));
    }

    // populate matrix
    while (p.hasNext())
    {
        MatrixElement e(p.getNext());
        this->operator()(e.row(), e.col()) = e.value();
    }
}

ColMatrix::ColMatrix(FileParser &p, bool partitionRows, std::vector<unsigned> whichIndices)
{
    mNumRows = whichIndices.size();
    mNumCols = partitionRows ? p.nCol() : p.nRow();

    // allocate matrix
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(Vector(mNumRows));
    }

    // fill in matrix
    fill(*this, p, partitionRows, whichIndices);
}

ColMatrix ColMatrix::operator*(float val) const
{
    ColMatrix mat(mNumRows, mNumCols);
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mat.getCol(j) = mCols[j] * val;
    }
    return mat;
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

ColMatrix& ColMatrix::operator=(const RowMatrix &mat)
{
    copyMatrix(*this, mat);
    return *this;
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
