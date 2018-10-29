#ifndef __COGAPS_VECTOR_H__
#define __COGAPS_VECTOR_H__

#include "../utils/Archive.h"

#include <boost/align/aligned_allocator.hpp>
#include <vector>

// need to align data for SIMD
namespace bal = boost::alignment;
typedef std::vector<float, bal::aligned_allocator<float,32> > aligned_vector;

// no iterator access, only random access
class Vector
{
public:

    explicit Vector(unsigned sz);
    explicit Vector(const std::vector<float> &v);

    float operator[](unsigned i) const;
    float& operator[](unsigned i);

    const float* ptr() const;
    float* ptr();

    unsigned size() const;
    void pad(float val);

    void operator+=(const Vector &v);

    void operator*=(float f);
    void operator/=(float f);
    
    friend Archive& operator<<(Archive &ar, const Vector &vec);
    friend Archive& operator>>(Archive &ar, Vector &vec);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    aligned_vector mData;
    unsigned mSize;
};

#endif