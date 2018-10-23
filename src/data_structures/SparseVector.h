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

    template <unsigned N>
    friend class TemplatedSparseIterator;
    
    friend class SparseMatrix; // for inserting values

    explicit SparseVector(unsigned size);
    explicit SparseVector(const std::vector<float> &v);
    explicit SparseVector(const Vector &v);

    unsigned size() const;

    Vector getDense() const;
    
    float at(unsigned n) const;
    unsigned nElements() const;

    friend Archive& operator<<(Archive &ar, const SparseVector &vec);
    friend Archive& operator>>(Archive &ar, SparseVector &vec);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif
    
    unsigned mSize;
    std::vector<uint64_t> mIndexBitFlags;
    std::vector<float> mData;

    void insert(unsigned i, float v);
};

#endif // __COGAPS_SPARSE_VECTOR_H__