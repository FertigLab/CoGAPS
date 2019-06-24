#ifndef __COGAPS_SPARSE_VECTOR_H__
#define __COGAPS_SPARSE_VECTOR_H__

#include <stdint.h>
#include <vector>

class Archive;
class Vector;
class SparseMatrix;

class SparseVector
{
public:
    explicit SparseVector(unsigned size);
    explicit SparseVector(const std::vector<float> &v);
    explicit SparseVector(const Vector &v);
    unsigned size() const;
    Vector getDense() const;
    const std::vector<float>& getData() const { return mData; }
    const std::vector<uint64_t>& getBitFlags() const { return mIndexBitFlags; }
    float at(unsigned n) const;
    float getIthElement(unsigned n) const;
    unsigned nElements() const;
    friend Archive& operator<<(Archive &ar, const SparseVector &vec);
    friend Archive& operator>>(Archive &ar, SparseVector &vec);
private:
    template <unsigned N>
    friend class SparseIterator;
    friend class SparseMatrix; // for inserting values
    unsigned mSize;
    std::vector<uint64_t> mIndexBitFlags;
    std::vector<float> mData;
    void insert(unsigned i, float v);
};

#endif // __COGAPS_SPARSE_VECTOR_H__