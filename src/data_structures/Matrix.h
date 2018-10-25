#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include "../utils/Archive.h"
#include "Vector.h"

#include <vector>

class Matrix
{
public:

    Matrix();
    Matrix(unsigned nrow, unsigned ncol);
    Matrix(const Matrix &mat, bool transpose, bool partitionRows,
        std::vector<unsigned> indices);
    Matrix(const std::string &path, bool transpose, bool partitionRows,
        std::vector<unsigned> indices);

    unsigned nRow() const;
    unsigned nCol() const;

    void pad(float val);

    float operator()(unsigned i, unsigned j) const;
    float& operator()(unsigned i, unsigned j);

    Vector& getCol(unsigned col);
    const Vector& getCol(unsigned col) const;
    
    bool empty() const;

    friend Archive& operator<<(Archive &ar, const Matrix &vec);
    friend Archive& operator>>(Archive &ar, Matrix &vec);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    std::vector<Vector> mCols;
    unsigned mNumRows;
    unsigned mNumCols;
};

#endif // __COGAPS_MATRIX_H__