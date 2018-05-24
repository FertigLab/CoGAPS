#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include "../Archive.h"
#include "Vector.h"

#include <Rcpp.h>
#include <vector>

struct MatrixChange
{
    char label;
    unsigned nChanges;

    unsigned row1;
    unsigned col1;
    float delta1;
    
    unsigned row2;
    unsigned col2;
    float delta2;

    MatrixChange(char l, unsigned r, unsigned c, float d)
        : label(l), nChanges(1), row1(r), col1(c), delta1(d), row2(0),
        col2(0), delta2(0.f)
    {}

    MatrixChange(char l, unsigned r1, unsigned c1, float d1, unsigned r2,
    unsigned c2, float d2)
        : label(l), nChanges(2), row1(r1), col1(c1), delta1(d1), row2(r2),
        col2(c2), delta2(d2)
    {}
};


class ColMatrix;

class RowMatrix
{
private:
    
    std::vector<Vector> mRows;
    unsigned mNumRows, mNumCols;

public:

    RowMatrix(unsigned nrow, unsigned ncol);
    RowMatrix(const Rcpp::NumericMatrix &rmat);
    RowMatrix(const std::string &path);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    float& operator()(unsigned r, unsigned c) {return mRows[r][c];}
    float operator()(unsigned r, unsigned c) const {return mRows[r][c];}

    Vector& getRow(unsigned row) {return mRows[row];}
    const Vector& getRow(unsigned row) const {return mRows[row];}

    const float* rowPtr(unsigned row) const {return mRows[row].ptr();}
    float* rowPtr(unsigned row) {return mRows[row].ptr();}

    RowMatrix operator/(float val) const;

    void operator=(const RowMatrix &mat);
    void operator=(const ColMatrix &mat);

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
    ColMatrix(const Rcpp::NumericMatrix &rmat);
    ColMatrix(const std::string &path);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    float& operator()(unsigned r, unsigned c) {return mCols[c][r];}
    float operator()(unsigned r, unsigned c) const {return mCols[c][r];}

    Vector& getCol(unsigned col) {return mCols[col];}
    const Vector& getCol(unsigned col) const {return mCols[col];}

    const float* colPtr(unsigned col) const {return mCols[col].ptr();}
    float* colPtr(unsigned col) {return mCols[col].ptr();}

    ColMatrix operator/(float val) const;

    void operator=(const ColMatrix &mat);
    void operator=(const RowMatrix &mat);

    Rcpp::NumericMatrix rMatrix() const;

    friend Archive& operator<<(Archive &ar, ColMatrix &mat);
    friend Archive& operator>>(Archive &ar, ColMatrix &mat);
};

#endif