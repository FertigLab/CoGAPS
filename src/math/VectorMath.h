#ifndef __COGAPS_VECTOR_MATH_H__
#define __COGAPS_VECTOR_MATH_H__

#include "../data_structures/Vector.h"
#include "../data_structures/HybridVector.h"
#include "../data_structures/SparseVector.h"
#include "../utils/GapsAssert.h"
#include "SIMD.h"

namespace gaps
{
    float min(const Vector &v);
    float min(const HybridVector &v);
    float min(const SparseVector &v);
    float max(const Vector &v);
    float max(const HybridVector &v);
    float max(const SparseVector &v);
    unsigned whichMax(const Vector &v);
    float sum(const Vector &v);
    float sum(const HybridVector &v);
    float sum(const SparseVector &v);
    bool isVectorZero(const Vector &v);
    bool isVectorZero(const HybridVector &v);
    Vector elementSq(Vector v);
    Vector pmax(Vector v, float p);
    template <class VectorType>
    float dot(const VectorType &a, const VectorType &b);
    template <class VectorType>
    float dot_diff(const VectorType &a, const VectorType &b, const VectorType &c);
} // namespace gaps

Vector operator*(Vector v, float f);
Vector operator/(Vector v, float f);
Vector operator*(const HybridVector &hv, float f);
Vector operator/(const HybridVector &hv, float f);

// this function is frequently called small vectors and represents a significant bottleneck,
// this code falls through the switch statement to avoid the overhead of branching,
// supported vectors have to be padded with 0 so that the length is a multiple of SIMD_INC
template <class VectorType>
float gaps::dot(const VectorType &a, const VectorType &b)
{
    GAPS_ASSERT(a.size() == b.size());
    GAPS_ASSERT((a.size() % SIMD_INC) == 0);
    const float *v1 = a.ptr();
    const float *v2 = b.ptr();
    const unsigned size = a.size();
    unsigned nChunks = 1 + (size - 1) / SIMD_INC;
    gaps_packed_t pDot(SET_SCALAR(0.f));
    switch (nChunks)
    {
        case 25:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 24 * SIMD_INC), LOAD_PACKED(v2 + 24 * SIMD_INC)));
            // fall through
        case 24:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 23 * SIMD_INC), LOAD_PACKED(v2 + 23 * SIMD_INC)));
            // fall through
        case 23:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 22 * SIMD_INC), LOAD_PACKED(v2 + 22 * SIMD_INC)));
            // fall through
        case 22:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 21 * SIMD_INC), LOAD_PACKED(v2 + 21 * SIMD_INC)));
            // fall through
        case 21:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 20 * SIMD_INC), LOAD_PACKED(v2 + 20 * SIMD_INC)));
            // fall through
        case 20:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 19 * SIMD_INC), LOAD_PACKED(v2 + 19 * SIMD_INC)));
            // fall through
        case 19:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 18 * SIMD_INC), LOAD_PACKED(v2 + 18 * SIMD_INC)));
            // fall through
        case 18:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 17 * SIMD_INC), LOAD_PACKED(v2 + 17 * SIMD_INC)));
            // fall through
        case 17:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 16 * SIMD_INC), LOAD_PACKED(v2 + 16 * SIMD_INC)));
            // fall through
        case 16:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 15 * SIMD_INC), LOAD_PACKED(v2 + 15 * SIMD_INC)));
            // fall through
        case 15:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 14 * SIMD_INC), LOAD_PACKED(v2 + 14 * SIMD_INC)));
            // fall through
        case 14:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 13 * SIMD_INC), LOAD_PACKED(v2 + 13 * SIMD_INC)));
            // fall through
        case 13:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 12 * SIMD_INC), LOAD_PACKED(v2 + 12 * SIMD_INC)));
            // fall through
        case 12:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 11 * SIMD_INC), LOAD_PACKED(v2 + 11 * SIMD_INC)));
            // fall through
        case 11:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 10 * SIMD_INC), LOAD_PACKED(v2 + 10 * SIMD_INC)));
            // fall through
        case 10:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 9 * SIMD_INC), LOAD_PACKED(v2 + 9 * SIMD_INC)));
            // fall through
        case 9:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 8 * SIMD_INC), LOAD_PACKED(v2 + 8 * SIMD_INC)));
            // fall through
        case 8:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 7 * SIMD_INC), LOAD_PACKED(v2 + 7 * SIMD_INC)));
            // fall through
        case 7:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 6 * SIMD_INC), LOAD_PACKED(v2 + 6 * SIMD_INC)));
            // fall through
        case 6:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 5 * SIMD_INC), LOAD_PACKED(v2 + 5 * SIMD_INC)));
            // fall through
        case 5:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 4 * SIMD_INC), LOAD_PACKED(v2 + 4 * SIMD_INC)));
            // fall through
        case 4:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 3 * SIMD_INC), LOAD_PACKED(v2 + 3 * SIMD_INC)));
            // fall through
        case 3:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + 2 * SIMD_INC), LOAD_PACKED(v2 + 2 * SIMD_INC)));
            // fall through
        case 2:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + SIMD_INC), LOAD_PACKED(v2 + SIMD_INC)));
            // fall through
        case 1:
            pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1), LOAD_PACKED(v2)));
            break;    
        default:
            for (unsigned i = 0; i < size; i += SIMD_INC)
            {
                pDot = ADD_PACKED(pDot, MUL_PACKED(LOAD_PACKED(v1 + i), LOAD_PACKED(v2 + i)));
            }
            break;
    }
    return gaps::simd::getScalar(pDot);
}

template <class VectorType>
float gaps::dot_diff(const VectorType &a, const VectorType &b, const VectorType &c)
{
    GAPS_ASSERT(a.size() == b.size());
    GAPS_ASSERT(a.size() == c.size());
    const float *v1 = a.ptr();
    const float *v2 = b.ptr();
    const float *v3 = c.ptr();
    const unsigned size = a.size();
    gaps::simd::PackedFloat packedDot(0.f), p1, p2, p3;
    for (unsigned i = 0; i < size; i += SIMD_INC)
    {
        p1.load(v1 + i);
        p2.load(v2 + i);
        p3.load(v3 + i);
        packedDot += p1 * (p2 - p3);
    }
    return packedDot.scalar();
}

#endif // __COGAPS_VECTOR_MATH_H__