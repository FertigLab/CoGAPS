#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include "Vector.h"

#include <string>
#include <vector>

class Archive;

class Matrix
{
public:

    Matrix();
    Matrix(unsigned nrow, unsigned ncol);
    Matrix(const Matrix &mat, bool genesInCols, bool subsetGenes,
        std::vector<unsigned> indices);
    Matrix(const std::string &path, bool genesInCols, bool subsetGenes,
        std::vector<unsigned> indices);
    unsigned nRow() const;
    unsigned nCol() const;
    void pad(float val);
    float operator()(unsigned i, unsigned j) const;
    float& operator()(unsigned i, unsigned j);
    Vector& getCol(unsigned col);
    const Vector& getCol(unsigned col) const;
    bool empty() const;
    Matrix getMatrix() const;
    friend Archive& operator<<(Archive &ar, const Matrix &mat);
    friend Archive& operator>>(Archive &ar, Matrix &mat);
private:
    std::vector<Vector> mCols;
    unsigned mNumRows;
    unsigned mNumCols;
};

#endif // __COGAPS_MATRIX_H__