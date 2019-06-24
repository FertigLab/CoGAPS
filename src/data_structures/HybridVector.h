#ifndef __COGAPS_HYBRID_VECTOR_H__
#define __COGAPS_HYBRID_VECTOR_H__

#include "../utils/Archive.h"

#include "boost/align/aligned_allocator.hpp"
#include <vector>

// need to align data for SIMD
namespace bal = boost::alignment;
typedef std::vector<float, bal::aligned_allocator<float,32> > aligned_vector;

class SparseIteratorTwo;
class SparseIteratorThree;

// stored as a dense vector (efficient setting of values) but maintains
// index bit flags of non-zeros so it can be used with SparseIterator
class HybridVector
{
public:

    friend class SparseIteratorTwo;
    friend class SparseIteratorThree;
    
    template <unsigned N>
    friend class SparseIterator;

    explicit HybridVector(unsigned sz);
    explicit HybridVector(const std::vector<float> &v);

    bool empty() const;
    unsigned size() const;

    bool add(unsigned i, float v); // true if zeros out data
    bool set(unsigned i, float v); // true if zeros out data
    float operator[](unsigned i) const;

    const float* ptr() const;

    const std::vector<uint64_t>& getBitFlags() const { return mIndexBitFlags; }

    friend Archive& operator<<(Archive &ar, const HybridVector &vec);
    friend Archive& operator>>(Archive &ar, HybridVector &vec);

private:

    std::vector<uint64_t> mIndexBitFlags;
    aligned_vector mData;
    unsigned mSize;
};

#endif // __COGAPS_HYBRID_VECTOR_H__