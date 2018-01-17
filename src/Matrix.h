#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include "Archive.h"

#include <Rcpp.h>
#include <vector>

// use CRTP for different matrix types?

// temporary: used for testing performance of float vs double
typedef float matrix_data_t;

struct MatrixChange
{
    char label;
    unsigned nChanges;

    unsigned row1;
    unsigned col1;
    matrix_data_t delta1;
    
    unsigned row2;
    unsigned col2;
    matrix_data_t delta2;

    MatrixChange(char l, unsigned r, unsigned c, float d)
        : label(l), nChanges(1), row1(r), col1(c), delta1(d), row2(0),
        col2(0), delta2(0.0)
    {}

    MatrixChange(char l, unsigned r1, unsigned c1, float d1, unsigned r2,
    unsigned c2, float d2)
        : label(l), nChanges(2), row1(r1), col1(c1), delta1(d1), row2(r2),
        col2(c2), delta2(d2)
    {}
};

class Vector
{
private:

    std::vector<matrix_data_t> mValues;

public:

    Vector(unsigned size) : mValues(std::vector<matrix_data_t>(size, 0.0)) {}
    Vector(const std::vector<matrix_data_t>& v) : mValues(v) {}
    Vector(const Vector &vec) : mValues(vec.mValues) {}

    matrix_data_t& operator[](unsigned i) {return mValues[i];}
    matrix_data_t operator[](unsigned i) const {return mValues[i];}
    unsigned size() const {return mValues.size();}

    Rcpp::NumericVector rVec() const {return Rcpp::wrap(mValues);}
    void concat(const Vector& vec);
    void operator+=(const Vector &vec);

    friend void operator<<(Archive &ar, Vector &vec);
    friend void operator>>(Archive &ar, Vector &vec);
};

class RowMatrix
{
private:
    
    std::vector<Vector> mRows;
    unsigned mNumRows, mNumCols;

public:

    RowMatrix(unsigned nrow, unsigned ncol);
    RowMatrix(const Rcpp::NumericMatrix& rmat);
    RowMatrix(const RowMatrix &mat);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    matrix_data_t& operator()(unsigned r, unsigned c) {return mRows[r][c];}
    matrix_data_t operator()(unsigned r, unsigned c) const {return mRows[r][c];}

    Vector& getRow(unsigned row) {return mRows[row];}
    const Vector& getRow(unsigned row) const {return mRows[row];}

    void update(const MatrixChange &change);
    Rcpp::NumericMatrix rMatrix() const;

    friend void operator<<(Archive &ar, RowMatrix &mat);
    friend void operator>>(Archive &ar, RowMatrix &mat);
};

class ColMatrix
{
private:

    std::vector<Vector> mCols;
    unsigned mNumRows, mNumCols;

public:

    ColMatrix(unsigned nrow, unsigned ncol);
    ColMatrix(const Rcpp::NumericMatrix& rmat);
    ColMatrix(const ColMatrix &mat);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    matrix_data_t& operator()(unsigned r, unsigned c) {return mCols[c][r];}
    matrix_data_t operator()(unsigned r, unsigned c) const {return mCols[c][r];}

    Vector& getCol(unsigned col) {return mCols[col];}
    const Vector& getCol(unsigned col) const {return mCols[col];}

    void update(const MatrixChange &change);
    Rcpp::NumericMatrix rMatrix() const;

    friend void operator<<(Archive &ar, ColMatrix &mat);
    friend void operator>>(Archive &ar, ColMatrix &mat);
};

// gain performance at the expense of memory
class TwoWayMatrix
{
private:

    RowMatrix mRowMatrix;
    ColMatrix mColMatrix;

public:

    TwoWayMatrix(unsigned nrow, unsigned ncol)
        : mRowMatrix(nrow, ncol), mColMatrix(nrow, ncol)
    {}

    TwoWayMatrix(const Rcpp::NumericMatrix& rmat)
        : mRowMatrix(rmat), mColMatrix(rmat)
    {}

    unsigned nRow() const {return mRowMatrix.nRow();}
    unsigned nCol() const {return mRowMatrix.nCol();}
    
    const Vector& getRow(unsigned row) const {return mRowMatrix.getRow(row);}
    const Vector& getCol(unsigned col) const {return mColMatrix.getCol(col);}

    void set(unsigned row, unsigned col, float value)
    {
        mRowMatrix(row,col) = value;
        mColMatrix(row,col) = value;
    }
    
    Rcpp::NumericMatrix rMatrix() const
    {
        return mRowMatrix.rMatrix();
    }

    friend void operator<<(Archive &ar, TwoWayMatrix &mat);
    friend void operator>>(Archive &ar, TwoWayMatrix &mat);
};

#endif