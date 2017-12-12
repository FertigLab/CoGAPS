#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include <Rcpp.h>
#include <vector>

// temporary: used for testing performance of float vs double
typedef double matrix_data_t;

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

    MatrixChange(char l, unsigned r, unsigned c, double d)
        : label(l), nChanges(1), row1(r), col1(c), delta1(d)
    {}

    MatrixChange(char l, unsigned r1, unsigned c1, double d1, unsigned r2,
    unsigned c2, double d2)
        : label(l), nChanges(2), row1(r1), col1(c1), delta1(d1), row2(r2),
        col2(c2), delta2(d2)
    {}
};

class Vector;
class RowMatrix;
class ColMatrix;
class TwoWayMatrix;

// no polymorphism to prevent virtual function overhead, not really
// needed anyways since few functions are used on all types of matrices

class Vector
{
private:

    std::vector<matrix_data_t> mValues;

public:

    Vector(unsigned size) : mValues(std::vector<matrix_data_t>(size, 0.0)) {}
    Vector(const std::vector<matrix_data_t>& v) : mValues(v) {}

    matrix_data_t& operator()(unsigned i) {return mValues[i];}
    matrix_data_t operator()(unsigned i) const {return mValues[i];}
    unsigned size() const {return mValues.size();}

    Rcpp::NumericVector rVec() const;
    void concat(const Vector& vec);

    void operator+=(const Vector &vec);
};

class RowMatrix
{
private:

    std::vector<Vector> mRows;
    unsigned mNumRows, mNumCols;

public:

    RowMatrix(unsigned nrow, unsigned ncol);
    RowMatrix(const Rcpp::NumericMatrix& rmat);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    matrix_data_t& operator()(unsigned r, unsigned c) {return mRows[r](c);}
    matrix_data_t operator()(unsigned r, unsigned c) const {return mRows[r](c);}

    Vector& getRow(unsigned row);
    const Vector& getRow(unsigned row) const;

    void update(const MatrixChange &change);
    Rcpp::NumericMatrix rMatrix() const;
};

class ColMatrix
{
private:

    std::vector<Vector> mCols;
    unsigned mNumRows, mNumCols;

public:

    ColMatrix(unsigned nrow, unsigned ncol);
    ColMatrix(const Rcpp::NumericMatrix& rmat);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    matrix_data_t& operator()(unsigned r, unsigned c) {return mCols[c](r);}
    matrix_data_t operator()(unsigned r, unsigned c) const {return mCols[c](r);}

    Vector& getCol(unsigned col);
    const Vector& getCol(unsigned col) const;

    void update(const MatrixChange &change);
    Rcpp::NumericMatrix rMatrix() const;
};

// gain performance at the cost of memory
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

    matrix_data_t operator()(unsigned r, unsigned c) const
    {
        return mRowMatrix(r,c);
    }

    const Vector& getRow(unsigned row) const {return mRowMatrix.getRow(row);}
    const Vector& getCol(unsigned col) const {return mColMatrix.getCol(col);}

    void set(unsigned row, unsigned col, double value)
    {
        mRowMatrix(row, col) = value;
        mColMatrix(row, col) = value;
    }
    
    Rcpp::NumericMatrix rMatrix() const
    {
        return mRowMatrix.rMatrix();
    }
};

#endif