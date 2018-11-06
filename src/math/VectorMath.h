#ifndef __COGAPS_VECTOR_MATH_H__
#define __COGAPS_VECTOR_MATH_H__

#include "../data_structures/Vector.h"
#include "../data_structures/HybridVector.h"

namespace gaps
{
    float min(const Vector &v);
    float min(const HybridVector &v);

    float max(const Vector &v);
    float max(const HybridVector &v);

    bool isVectorZero(const Vector &v);
    bool isVectorZero(const HybridVector &v);

    float dot(const Vector &v1, const Vector &v2);
    float dot(const HybridVector &v1, const HybridVector &v2);

    float sum(const Vector &v);
    float sum(const HybridVector &v);

    Vector elementSq(Vector v);
    Vector pmax(Vector v, float p);
} // namespace gaps

Vector operator*(Vector v, float f);
Vector operator/(Vector v, float f);

Vector operator*(const HybridVector &hv, float f);
Vector operator/(const HybridVector &hv, float f);

#endif