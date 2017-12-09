#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include <Rcpp.h>
#include <vector>

// temporary: used for testing performance of float vs double
typedef double matrix_data_t;

struct MatrixChange
{
    char type;
    unsigned nChanges;

    unsigned row1;
    unsigned col1;
    matrix_data_t delta1;
    
    unsigned row2;
    unsigned col2;
    matrix_data_t delta2;

    MatrixChange(char t, unsigned r, unsigned c, double d)
        : type(t), nChanges(1), row1(r), col1(c), delta1(d)
    {}

    MatrixChange(char t, unsigned r1, unsigned c1, double d1, unsigned r2, unsigned c2, double d2)
        : type(t), nChanges(2), row1(r1), col1(c1), delta1(d1), row2(r2), col2(c2), delta2(d2)
    {}
};

class Vector;
class RowMatrix;
class ColMatrix;
typedef RowMatrix Matrix; // default matrix type

class Vector
{
private:

    std::vector<matrix_data_t> mValues;

public:

    Vector(unsigned size) : mValues(std::vector<matrix_data_t>(size, 0)) {}
    Vector(const std::vector<matrix_data_t>& v) : mValues(v) {}

    Rcpp::NumericVector rVector() const {return Rcpp::wrap(mValues);}

    matrix_data_t& operator()(unsigned int i) {return mValues[i];}
    matrix_data_t operator()(unsigned int i) const {return mValues[i];}
    unsigned size() const {return mValues.size();}

    matrix_data_t dotProduct(const Vector &vec) const;
    matrix_data_t sum() const;

    typedef std::vector<matrix_data_t>::iterator iterator;
    iterator begin() {return mValues.begin();}
    iterator end() {return mValues.end();}
};

class RowMatrix
{
private:

    std::vector<Vector> mRows;
    unsigned int mNumRows, mNumCols;

public:

    RowMatrix(unsigned int nrow, unsigned int ncol);
    RowMatrix(const std::vector< std::vector<matrix_data_t> > &mat);
    RowMatrix(const Rcpp::NumericMatrix& rmat);

    unsigned int nRow() const {return mNumRows;}
    unsigned int nCol() const {return mNumCols;}

    matrix_data_t& operator()(unsigned int r, unsigned int c) {return mRows[r](c);}
    matrix_data_t operator()(unsigned int r, unsigned int c) const {return mRows[r](c);}

    Vector& getRow(unsigned int row);
    const Vector& getRow(unsigned int row) const;

    void elemUpdate(const MatrixChange &change);
};

class ColMatrix
{
private:

    std::vector<Vector> mCols;
    unsigned int mNumRows, mNumCols;

public:

    ColMatrix(unsigned int nrow, unsigned int ncol);
    ColMatrix(const Rcpp::NumericMatrix& rmat);

    unsigned int nRow() const {return mNumRows;}
    unsigned int nCol() const {return mNumCols;}

    matrix_data_t& operator()(unsigned int r, unsigned int c) {return mCols[c](r);}
    matrix_data_t operator()(unsigned int r, unsigned int c) const {return mCols[c](r);}

    Vector& getCol(unsigned int col);

    void elemUpdate(const MatrixChange &change);
};

#endif