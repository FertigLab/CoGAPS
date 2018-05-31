#ifndef __COGAPS_VECTOR_H__
#define __COGAPS_VECTOR_H__

#include "../Archive.h"

#include <boost/align/aligned_allocator.hpp>
#include <vector>

// need to align data for SIMD
namespace bal = boost::alignment;
typedef std::vector<float, bal::aligned_allocator<float,32> > aligned_vector;

class Vector
{
private:

    aligned_vector mValues;

public:

    explicit Vector(unsigned size) : mValues(aligned_vector(size, 0.f)) {}
    explicit Vector(const std::vector<float> &v);

    const float* ptr() const {return &mValues[0];}
    float* ptr() {return &mValues[0];}

    float& operator[](unsigned i) {return mValues[i];}
    float operator[](unsigned i) const {return mValues[i];}

    unsigned size() const {return mValues.size();}

    Rcpp::NumericVector rVec() const {return Rcpp::wrap(mValues);}
    void concat(const Vector& vec);
    void operator+=(const Vector &vec);

    Vector operator*(float val) const;
    Vector operator/(float val) const;

    friend Archive& operator<<(Archive &ar, Vector &vec);
    friend Archive& operator>>(Archive &ar, Vector &vec);
};

#endif