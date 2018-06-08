#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include "../Archive.h"
#include "../file_parser/MatrixElement.h"
#include "Vector.h"

#include <Rcpp.h>
#include <vector>
#include <algorithm>

// forward declarations
class RowMatrix;
class ColMatrix;

class RowMatrix
{
private:

    std::vector<Vector> mRows;
    unsigned mNumRows, mNumCols;

public:

    RowMatrix(unsigned nrow, unsigned ncol);
    explicit RowMatrix(const Rcpp::NumericMatrix &rmat);

    template <class Parser>
    RowMatrix(Parser &p, bool parseRows, std::vector<unsigned> whichIndices);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    float& operator()(unsigned r, unsigned c) {return mRows[r][c];}
    float operator()(unsigned r, unsigned c) const {return mRows[r][c];}

    Vector& getRow(unsigned row) {return mRows[row];}
    const Vector& getRow(unsigned row) const {return mRows[row];}

    const float* rowPtr(unsigned row) const {return mRows[row].ptr();}
    float* rowPtr(unsigned row) {return mRows[row].ptr();}

    RowMatrix operator/(float val) const;

    RowMatrix& operator=(const ColMatrix &mat);

    Rcpp::NumericMatrix rMatrix() const;

    friend Archive& operator<<(Archive &ar, RowMatrix &mat);
    friend Archive& operator>>(Archive &ar, RowMatrix &mat);
};

class ColMatrix
{
private:

    std::vector<Vector> mCols;
    unsigned mNumRows, mNumCols;

public:

    ColMatrix(unsigned nrow, unsigned ncol);
    explicit ColMatrix(const Rcpp::NumericMatrix &rmat);

    template <class Parser>
    ColMatrix(Parser &p, bool parseRows, std::vector<unsigned> whichIndices);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    float& operator()(unsigned r, unsigned c) {return mCols[c][r];}
    float operator()(unsigned r, unsigned c) const {return mCols[c][r];}

    Vector& getCol(unsigned col) {return mCols[col];}
    const Vector& getCol(unsigned col) const {return mCols[col];}

    const float* colPtr(unsigned col) const {return mCols[col].ptr();}
    float* colPtr(unsigned col) {return mCols[col].ptr();}

    ColMatrix operator/(float val) const;

    ColMatrix& operator=(const RowMatrix &mat);

    Rcpp::NumericMatrix rMatrix() const;

    friend Archive& operator<<(Archive &ar, ColMatrix &mat);
    friend Archive& operator>>(Archive &ar, ColMatrix &mat);
};


// if partitionRows is false, partition columns instead
// rows of matrix should be partition dimension, i.e. need to transpose
// is partitionRows is false
template <class Matrix, class Parser>
inline fill(Matrix &mat, Parser &p, bool partitionRows, std::vector<unsigned> whichIndices)
{
    // TODO implement
}

template <class Parser>
RowMatrix::RowMatrix(Parser &p, bool partitionRows, std::vector<unsigned> whichIndices)
: mNumRows(?), mNumCols(?)
{
    // allocate matrix
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
    }

    // fill in matrix
    fill(*this, p, partitionRows, whichIndices);
}

template <class Parser>
ColMatrix::ColMatrix(Parser &p, bool partitionRows, std::vector<unsigned> whichIndices)
: mNumRows(?), mNumCols(?)
{
    // allocate matrix
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(Vector(mNumRows));
    }

    // fill in matrix
    fill(*this, p, partitionRows, whichIndices);
}


//// BELOW CODE IS OUTDATED, REMOVE BEFORE MERGE

// construct RowMatrix from file
template <class Parser>
RowMatrix::RowMatrix(Parser &p, unsigned nrow, unsigned ncol) : mNumRows(nrow),
mNumCols(ncol)
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
        this->operator()(e.row, e.col) = e.value;
    }
}

// construct ColMatrix from file
template <class Parser>
ColMatrix::ColMatrix(Parser &p, unsigned nrow, unsigned ncol) : mNumRows(nrow),
mNumCols(ncol)
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
        this->operator()(e.row, e.col) = e.value;
    }
}

// This should construct a matrix, only reading the columns in "whichCols"
// from the file. The matrix should have "nrow" rows and "whichCols.size()"
// columns.
template <class Parser>
RowMatrix::RowMatrix(Parser &p, unsigned nrow, std::vector<unsigned> whichCols)
{
    // TODO implement
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(whichCols.size()));
    }

    while (p.hasNext())
    {
        MatrixElement e(p.getNext());
        std::vector<unsigned>::iterator newColsIndex = std::find(whichCols.begin(), whichCols.end(), e.col);
        if (newColsIndex != whichCols.end())
        {
            this->operator()(e.row, std::distance(whichCols.begin(), newColsIndex)) = e.value;
        }
    }
}

template <class Parser>
ColMatrix::ColMatrix(Parser &p, unsigned nrow, std::vector<unsigned> whichCols)
{
    // TODO implement
    for (unsigned j = 0; j < whichCols.size(); ++j)
    {
        mCols.push_back(Vector(mNumRows));
    }

    while(p.hasNext())
    {
        MatrixElement e(p.getNext());
        std::vector<unsigned>::iterator newColsIndex = std::find(whichCols.begin(), whichCols.end(), e.col);
        if (newColsIndex != whichCols.end())
        {
            this->operator()(e.row, std::distance(whichCols.begin(), newColsIndex)) = e.value;
        }
    }
}

#endif
