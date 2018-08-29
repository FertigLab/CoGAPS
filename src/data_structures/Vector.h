#ifndef __COGAPS_VECTOR_H__
#define __COGAPS_VECTOR_H__

#include "../utils/Archive.h"

#include <boost/align/aligned_allocator.hpp>
#include <vector>

// need to align data for SIMD
namespace bal = boost::alignment;
typedef std::vector<float, bal::aligned_allocator<float,32> > aligned_vector;

class Vector
{
public:

    explicit Vector(unsigned size);
    explicit Vector(const std::vector<float> &v);

    unsigned size() const;

    float* ptr();
    const float* ptr() const;

    float& operator[](unsigned i);
    float operator[](unsigned i) const;

    void operator+=(const Vector &vec);
    Vector operator-(Vector v) const;
    Vector operator*(float val) const;
    Vector operator/(float val) const;

    void operator*=(float val);
    void operator/=(float val);

    friend Archive& operator<<(Archive &ar, Vector &vec);
    friend Archive& operator>>(Archive &ar, Vector &vec);

private:

    aligned_vector mValues;
};

#endif // __COGAPS_VECTOR_H__