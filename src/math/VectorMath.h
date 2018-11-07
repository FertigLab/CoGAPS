#ifndef __COGAPS_VECTOR_MATH_H__
#define __COGAPS_VECTOR_MATH_H__

#include "../data_structures/Vector.h"
#include "../data_structures/HybridVector.h"
#include "SIMD.h"

namespace gaps
{
    float min(const Vector &v);
    float min(const HybridVector &v);
    float min(const SparseVector &v);

    float max(const Vector &v);
    float max(const HybridVector &v);
    float max(const SparseVector &v);

    float sum(const Vector &v);
    float sum(const HybridVector &v);
    float sum(const SparseVector &v);

    bool isVectorZero(const Vector &v);
    bool isVectorZero(const HybridVector &v);

    template <class VectorType>
    float dot(const VectorType &a, const VectorType &b);

    template <class VectorType>
    float dot_shifted(const VectorType &v1, const VectorType &v2, float c);

    Vector elementSq(Vector v);
    Vector pmax(Vector v, float p);
} // namespace gaps

Vector operator*(Vector v, float f);
Vector operator/(Vector v, float f);

Vector operator*(const HybridVector &hv, float f);
Vector operator/(const HybridVector &hv, float f);

template <class VectorType>
float gaps::dot(const VectorType &a, const VectorType &b)
{
    GAPS_ASSERT(a.size() == b.size());

    const float *v1 = a.ptr();
    const float *v2 = b.ptr();
    const unsigned size = a.size();

    gaps::simd::PackedFloat packedDot(0.f), p1, p2;
    for (gaps::simd::Index i(0); i < size; ++i)
    {
        p1.load(v1 + i);
        p2.load(v2 + i);
        packedDot += p1 * p2;
    }
    return packedDot.scalar();
}

template <class VectorType>
float gaps::dot_shifted(const VectorType &a, const VectorType &b, float c)
{
    GAPS_ASSERT(a.size() == b.size());

    const float *v1 = a.ptr();
    const float *v2 = b.ptr();
    const unsigned size = a.size();

    gaps::simd::PackedFloat packedDot(0.f), shift(c), p1, p2;
    for (gaps::simd::Index i(0); i < size; ++i)
    {
        p1.load(v1 + i);
        p2.load(v2 + i);
        packedDot += p1 * (p2 + shift);
    }
    return packedDot.scalar();
}

#endif