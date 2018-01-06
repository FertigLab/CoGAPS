#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include <Rcpp.h>
#include <vector>
//#include <boost/serialization/vector.hpp>

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
        : label(l), nChanges(1), row1(r), col1(c), delta1(d), row2(0),
        col2(0), delta2(0.0)
    {}

    MatrixChange(char l, unsigned r1, unsigned c1, double d1, unsigned r2,
    unsigned c2, double d2)
        : label(l), nChanges(2), row1(r1), col1(c1), delta1(d1), row2(r2),
        col2(c2), delta2(d2)
    {}
};

// no polymorphism to prevent virtual function overhead, not really
// needed anyways since few functions are used on all types of matrices

class Vector
{
private:

    std::vector<matrix_data_t> mValues;

    /*friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar)
    {
        ar & mValues;
    }*/

public:

    Vector(unsigned size) : mValues(std::vector<matrix_data_t>(size, 0.0)) {}
    Vector(const std::vector<matrix_data_t>& v) : mValues(v) {}
    Vector(const Vector &vec) : mValues(vec.mValues) {}

    matrix_data_t& operator()(unsigned i) {return mValues[i];}
    matrix_data_t operator()(unsigned i) const {return mValues[i];}
    matrix_data_t& operator[](unsigned i) {return mValues[i];}
    matrix_data_t operator[](unsigned i) const {return mValues[i];}
    unsigned size() const {return mValues.size();}

    Rcpp::NumericVector rVec() const {return Rcpp::wrap(mValues);}
    void concat(const Vector& vec);
    void operator+=(const Vector &vec);
};

class RowMatrix
{
private:
    
    std::vector<Vector> mRows;
    unsigned mNumRows, mNumCols;

    /*friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar)
    {
        ar & mRows;
        ar & mNumRows;
        ar & mNumCols;
    }*/

public:

    RowMatrix(unsigned nrow, unsigned ncol);
    RowMatrix(const Rcpp::NumericMatrix& rmat);
    RowMatrix(const RowMatrix &mat);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    matrix_data_t& operator()(unsigned r, unsigned c) {return mRows[r](c);}
    matrix_data_t operator()(unsigned r, unsigned c) const {return mRows[r](c);}

    Vector& getRow(unsigned row) {return mRows[row];}
    const Vector& getRow(unsigned row) const {return mRows[row];}

    void update(const MatrixChange &change);
    Rcpp::NumericMatrix rMatrix() const;
};

class ColMatrix
{
private:

    std::vector<Vector> mCols;
    unsigned mNumRows, mNumCols;

    /*friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar)
    {
        ar & mCols;
        ar & mNumRows;
        ar & mNumCols;
    }*/

public:

    ColMatrix(unsigned nrow, unsigned ncol);
    ColMatrix(const Rcpp::NumericMatrix& rmat);
    ColMatrix(const ColMatrix &mat);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    matrix_data_t& operator()(unsigned r, unsigned c) {return mCols[c](r);}
    matrix_data_t operator()(unsigned r, unsigned c) const {return mCols[c](r);}

    Vector& getCol(unsigned col) {return mCols[col];}
    const Vector& getCol(unsigned col) const {return mCols[col];}

    void update(const MatrixChange &change);
    Rcpp::NumericMatrix rMatrix() const;
};

// gain performance at the expense of memory
class TwoWayMatrix
{
private:

    RowMatrix mRowMatrix;
    ColMatrix mColMatrix;

    /*friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar)
    {
        ar & mRowMatrix;
        ar & mColMatrix;
    }*/

public:

    TwoWayMatrix(unsigned nrow, unsigned ncol)
        : mRowMatrix(nrow, ncol), mColMatrix(nrow, ncol)
    {}

    TwoWayMatrix(const Rcpp::NumericMatrix& rmat)
        : mRowMatrix(rmat), mColMatrix(rmat)
    {}

    unsigned nRow() const {return mRowMatrix.nRow();}
    unsigned nCol() const {return mRowMatrix.nCol();}
    
    // TODO remove since accessing this way defeats the purpose
    /*matrix_data_t operator()(unsigned r, unsigned c) const
    {
        return mRowMatrix(r,c);
    }*/

    const Vector& getRow(unsigned row) const {return mRowMatrix.getRow(row);}
    const Vector& getCol(unsigned col) const {return mColMatrix.getCol(col);}

    void set(unsigned row, unsigned col, double value)
    {
        mRowMatrix(row,col) = value;
        mColMatrix(row,col) = value;
    }
    
    Rcpp::NumericMatrix rMatrix() const
    {
        return mRowMatrix.rMatrix();
    }
};

#endif