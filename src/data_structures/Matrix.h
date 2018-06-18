#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include "../Archive.h"
#include "../file_parser/FileParser.h"
#include "Vector.h"

#include <algorithm>
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

    RowMatrix(const std::string &p);
    RowMatrix(const ColMatrix &mat);
    //RowMatrix(const std::string &p, bool parseRows, std::vector<unsigned> whichIndices);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    float& operator()(unsigned r, unsigned c) {return mRows[r][c];}
    float operator()(unsigned r, unsigned c) const {return mRows[r][c];}

    Vector& getRow(unsigned row) {return mRows[row];}
    const Vector& getRow(unsigned row) const {return mRows[row];}

    float* rowPtr(unsigned row) {return mRows[row].ptr();}
    const float* rowPtr(unsigned row) const {return mRows[row].ptr();}

    RowMatrix operator*(float val) const;
    RowMatrix operator/(float val) const;
    RowMatrix& operator=(const ColMatrix &mat);
    RowMatrix pmax(float scale) const;

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
    ColMatrix(const RowMatrix &mat);

    ColMatrix(const std::string &p);
    //ColMatrix(const std::string &p, bool parseRows, std::vector<unsigned> whichIndices);

    unsigned nRow() const {return mNumRows;}
    unsigned nCol() const {return mNumCols;}

    float& operator()(unsigned r, unsigned c) {return mCols[c][r];}
    float operator()(unsigned r, unsigned c) const {return mCols[c][r];}

    Vector& getCol(unsigned col) {return mCols[col];}
    const Vector& getCol(unsigned col) const {return mCols[col];}

    float* colPtr(unsigned col) {return mCols[col].ptr();}
    const float* colPtr(unsigned col) const {return mCols[col].ptr();}

    ColMatrix operator*(float val) const;
    ColMatrix operator/(float val) const;
    ColMatrix& operator=(const RowMatrix &mat);
    ColMatrix pmax(float scale) const;

    friend Archive& operator<<(Archive &ar, ColMatrix &mat);
    friend Archive& operator>>(Archive &ar, ColMatrix &mat);
};

#endif
