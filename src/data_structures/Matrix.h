#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include "../utils/Archive.h"
#include "Vector.h"

#include <vector>

class ColMatrix;
typedef ColMatrix Matrix; // when we don't care about row/col major order

class ColMatrix
{
public:

    // empty constructor
    ColMatrix(unsigned nrow, unsigned ncol);

    // construct with data
    ColMatrix(const ColMatrix &mat, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);
    ColMatrix(const std::string &path, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);

    unsigned nRow() const;
    unsigned nCol() const;

    ColMatrix operator*(float val) const;
    ColMatrix operator/(float val) const;

    // for single element access - do not loop over elements with this
    float& operator()(unsigned r, unsigned c);
    float operator()(unsigned r, unsigned c) const;

    // for convenience when doing non-performance critical math
    Vector& getCol(unsigned col);
    const Vector& getCol(unsigned col) const;

    // raw pointers are used for looping through a column with SIMD
    float* colPtr(unsigned col);
    const float* colPtr(unsigned col) const;

    friend Archive& operator<<(Archive &ar, ColMatrix &mat);
    friend Archive& operator>>(Archive &ar, ColMatrix &mat);

private:

    std::vector<Vector> mData;
    unsigned mNumRows, mNumCols;

    void allocate();
};

#endif // __COGAPS_MATRIX_H__