#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include "../Archive.h"
#include "../file_parser/MatrixElement.h"
#include "Vector.h"

#include <Rcpp.h>
#include <vector>

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
    RowMatrix(Parser &p, unsigned nrow, unsigned ncol);

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
    ColMatrix(Parser &p, unsigned nrow, unsigned ncol);

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

#endif