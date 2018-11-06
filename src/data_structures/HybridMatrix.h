#ifndef __COGAPS_HYBRID_MATRIX_H__
#define __COGAPS_HYBRID_MATRIX_H__

#include "HybridVector.h"
#include "Matrix.h"
#include "Vector.h"

#include "../utils/GapsAssert.h"

#include <vector>

// Stores data in row major order with dense vectors, and column major order
// with hybrid vectors - which provide a sparse vector interface but have
// a dense vector implementation. This data structure provides row-wise
// and column-wise dense vectors, while also allowing the columns to 
// be used in a SparseIterator with an actual SparseVector. Importantly,
// we can update values in O(1).
class HybridMatrix
{
public:

    HybridMatrix(unsigned nrow, unsigned ncol);
    
    unsigned nRow() const;
    unsigned nCol() const;

    void add(unsigned i, unsigned j, float v);
    float operator()(unsigned i, unsigned j) const;
    
    const Vector& getRow(unsigned n) const;
    const HybridVector& getCol(unsigned n) const;

    void operator=(const Matrix &mat);

    friend Archive& operator<<(Archive &ar, const HybridMatrix &vec);
    friend Archive& operator>>(Archive &ar, HybridMatrix &vec);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    std::vector<Vector> mRows;
    std::vector<HybridVector> mCols;
    unsigned mNumRows;
    unsigned mNumCols;
};

#endif // __COGAPS_HYBRID_MATRIX_H__