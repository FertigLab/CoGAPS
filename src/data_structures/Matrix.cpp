#include "Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/MatrixElement.h"
#include "../GapsAssert.h"

void RowMatrix::allocate()
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mData.push_back(Vector(mNumCols));
    }
}

float& RowMatrix::operator()(unsigned r, unsigned c)
{
    return mData[r][c];
}

float RowMatrix::operator()(unsigned r, unsigned c) const
{
    return mData[r][c];
}

Vector& getRow(unsigned row) {return mData[row];}
const Vector& getRow(unsigned row) const {return mData[row];}

float* rowPtr(unsigned row) {return mData[row].ptr();}
const float* rowPtr(unsigned row) const {return mData[row].ptr();}

void RowMatrix::allocate()
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
    }
}

float& operator()(unsigned r, unsigned c) {return mData[c][r];}
float operator()(unsigned r, unsigned c) const {return mData[c][r];}

Vector& getCol(unsigned col) {return mData[col];}
const Vector& getCol(unsigned col) const {return mData[col];}

float* colPtr(unsigned col) {return mData[col].ptr();}
const float* colPtr(unsigned col) const {return mData[col].ptr();}


template <class MatA, class MatB>
static void copyMatrix(MatA &dest, const MatB &source, bool transpose=false)
{
    GAPS_ASSERT(dest.nRow() == (transpose ? source.nCol() : source.nRow()));
    GAPS_ASSERT(dest.nCol() == (transpose ? source.nRow() : source.nCol()))
    for (unsigned i = 0; i < dest.nRow(); ++i)
    {
        for (unsigned j = 0; j < dest.nCol(); ++j)
        {
            dest(i,j) = transpose ? source(j,i) : source(i,j);
        }
    }
}

template <class Matrix>
static void createFromFile(Matrix &mat, const std::string &path, bool transpose)
{
    FileParser p(path);
    mNumRows = transpose ? p.nCol() : p.nRow();
    mNumCols = transpose ? p.nRow() : p.nCol();
    allocate();

    while (p.hasNext())
    {
        MatrixElement e(p.getNext());
        unsigned row = transpose ? e.col() : e.row();
        unsigned col = transpose ? e.row() : e.col();
        this->operator()(row, col) = e.value();
    }
}

// apply transpose first, then partition either rows or columns
static void createFromMatrix(const RowMatrix &mat, bool transpose, bool partitionRows,
const std::vector<unsigned> &indices)
{
    // create matrix
    mNumRows = partitionRows ? indices.size() : (transpose ? mat.nCol() : mat.nRow());
    mNumCols = partitionRows ? (transpose ? mat.nRow() : mat.nCol()) : indices.size();
    allocate();

    // fill matrix
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            unsigned dataIndex = transpose ? (partitionRows ? j : i) : (partitionRows ? i : j);
            std::vector<unsigned>::iterator pos = std::find(indices.begin(), indices.end(), dataIndex);
            if (pos != indices.end())
            {
                unsigned index = std::distance(indices.begin(), pos);
                unsigned row = partitionRows ? index : (transpose ? j : i);
                unsigned col = partitionRows ? (transpose ? i : j) : index;
                this->operator()(row,col) = mat(i,j)
            }
        }
    }
}

/****************************** ROW MATRIX *****************************/

void RowMatrix::allocate()
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
    }
}    

RowMatrix::RowMatrix(unsigned nrow, unsigned ncol)
: mNumRows(nrow), mNumCols(ncol)
{
    allocate();
}

RowMatrix::RowMatrix(const ColMatrix &mat)
: mNumRows(mat.nRow()), mNumCols(mat.nCol())
{
    allocate();
    copyMatrix(*this, mat);
}

RowMatrix::RowMatrix(const RowMatrix &mat, bool transpose)
    :
mNumRows(transpose ? mat.nCol() : mat.nRow()),
mNumCols(transpose ? mat.nRow() : mat.nCol())
{
    allocate();
    copyMatrix(*this, mat, true);
}

RowMatrix::RowMatrix(const std::string &path, bool transpose)
{
    createFromFile(path, transpose);
}

// apply transpose first, then partition either rows or columns
RowMatrix::RowMatrix(const RowMatrix &mat, bool transpose, bool partitionRows,
const std::vector<unsigned> &indices)
{
    // create matrix
    mNumRows = partitionRows ? indices.size() : (transpose ? mat.nCol() : mat.nRow());
    mNumCols = partitionRows ? (transpose ? mat.nRow() : mat.nCol()) : indices.size();
    allocate();

    // fill matrix
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            unsigned dataIndex = transpose ? (partitionRows ? j : i) : (partitionRows ? i : j);
            std::vector<unsigned>::iterator pos = std::find(indices.begin(), indices.end(), dataIndex);
            if (pos != indices.end())
            {
                unsigned index = std::distance(indices.begin(), pos);
                unsigned row = partitionRows ? index : (transpose ? j : i);
                unsigned col = partitionRows ? (transpose ? i : j) : index;
                this->operator()(row,col) = mat(i,j)
            }
        }
    }
}

// apply transpose first, then partition either rows or columns
RowMatrix::RowMatrix(const std::string &path, bool transpose, bool partitionRows,
const std::vector<unsigned> &indices)
{
    // create matrix
    FileParser p(path);
    mNumRows = partitionRows ? indices.size() : (transpose ? p.nCol() : p.nRow());
    mNumCols = partitionRows ? (transpose ? p.nRow() : p.nCol()) : indices.size();
    allocate();

    // fill
    unsigned indexLoc = transpose ? (partitionRows ? 1 : 0) : (partitionRows ? 0 : 1);
    unsigned rowLoc = transpose ? 1 : 0;
    unsigned colLoc = transpose ? 0 : 1;
    while (p.hasNext())
    {
        MatrixElement e(p.getNext());

        unsigned dataIndex = e[indexLoc];
        std::vector<unsigned>::iterator pos = std::find(indices.begin(), indices.end(), dataIndex);
        if (pos != indices.end())
        {
            unsigned index = std::distance(indices.begin(), pos);
            unsigned row = partitionRows ? index : e[rowLoc];
            unsigned col = partitionRows ? e[colLoc] : index;
            this->operator()(row, col) = e.value();
        }
    }
}


// if partitionRows is false, partition columns instead
// rows of matrix should be partition dimension, i.e. need to transpose
// is partitionRows is false
/*
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
            unsigned row = 
            mat.operator()(row, e[colSelect]) = e.value();
        }
    }
}
*/


/*
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
*/

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

RowMatrix RowMatrix::pmax(float scale, float max) const
{
    RowMatrix mat(mNumRows, mNumCols);
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            mat(i,j) = std::max(this->operator()(i,j) * scale, max);
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

ColMatrix::ColMatrix(const std::string &path)
{
    FileParser p(path);
    mNumRows = p.nRow();
    mNumCols = p.nCol();

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

/*
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
*/

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

ColMatrix ColMatrix::pmax(float scale, float max) const
{
    ColMatrix mat(mNumRows, mNumCols);
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            mat(i,j) = std::max(this->operator()(i,j) * scale, max);
        }
    }
    return mat;
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
