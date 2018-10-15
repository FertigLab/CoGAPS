#ifndef __COGAPS_SPARSE_VECTOR_H__
#define __COGAPS_SPARSE_VECTOR_H__

#include "Vector.h"
#include "../utils/Archive.h"

#include <vector>

class SparseMatrix;

class SparseVector
{
public:

    friend class SparseIterator;
    friend class SparseIteratorTwo;
    friend class SparseIteratorThree;
    friend class SparseMatrix; // for inserting values

    explicit SparseVector(unsigned size);
    explicit SparseVector(const std::vector<float> &v);
    explicit SparseVector(const Vector &v);

    unsigned size() const;

    Vector getDense() const;

    friend Archive& operator<<(Archive &ar, SparseVector &vec);
    friend Archive& operator>>(Archive &ar, SparseVector &vec);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif
    
    std::vector<uint64_t> mIndexBitFlags;
    std::vector<float> mData;
    unsigned mSize;

    void insert(unsigned i, float v);
};

#endif // __COGAPS_SPARSE_VECTOR_H__