#include "Math.h"
#include "VectorMath.h"
#include "SIMD.h"

static float dot_helper(const float *v1, const float *v2, unsigned size)
{
    gaps::simd::PackedFloat partialDot(0.f), p1, p2;
    gaps::simd::Index i(0);
    for (; i <= size - gaps::simd::Index::increment(); ++i)
    {
        p1.load(v1 + i);
        p2.load(v2 + i);
        partialDot += p1 * p2;
    }

    float dot = partialDot.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        dot += v1[j] * v2[j];
    }
    return dot;
}

float gaps::min(const Vector &v)
{
    float mn = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (v[i] < mn)
        {
            mn = v[i];
        }
    }
    return mn;
}

float gaps::min(const HybridVector &v)
{
    float mn = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (v[i] > mn)
        {
            mn = v[i];
        }
    }
    return mn;
}

float gaps::max(const Vector &v)
{
    float mx = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (v[i] > mx)
        {
            mx = v[i];
        }
    }
    return mx;
}

float gaps::max(const HybridVector &v)
{
    float mx = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (v[i] > mx)
        {
            mx = v[i];
        }
    }
    return mx;
}

bool gaps::isVectorZero(const Vector &v)
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (v[i] > 0.f)
        {
            return false;
        }
    }
    return true;
}

bool gaps::isVectorZero(const HybridVector &v)
{
    return v.empty();
}

float gaps::dot(const Vector &v1, const Vector &v2)
{
    GAPS_ASSERT(v1.size() == v2.size());

    return dot_helper(v1.ptr(), v2.ptr(), v1.size());
}

float gaps::dot(const HybridVector &v1, const HybridVector &v2)
{
    GAPS_ASSERT(v1.size() == v2.size());

    return dot_helper(v1.densePtr(), v2.densePtr(), v1.size());
}

float gaps::sum(const Vector &v)
{
    float sum = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        sum += v[i];
    }
    return sum;
}

float gaps::sum(const HybridVector &v)
{
    float sum = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        sum += v[i];
    }
    return sum;
}

Vector gaps::elementSq(Vector v)
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        v[i] *= v[i];
    }
    return v;
}

Vector operator*(Vector v, float f)
{
    v *= f;
    return v;
}

Vector operator/(Vector v, float f)
{
    v /= f;
    return v;
}

Vector operator*(const HybridVector &hv, float f)
{
    Vector v(hv.size());
    for (unsigned i = 0; i < hv.size(); ++i)
    {
        v[i] = hv[i] * f;
    }
    return v;
}

Vector operator/(const HybridVector &hv, float f)
{
    Vector v(hv.size());
    for (unsigned i = 0; i < hv.size(); ++i)
    {
        v[i] = hv[i] / f;
    }
    return v;
}

Vector gaps::pmax(Vector v, float p)
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        v[i] = gaps::max(v[i] * p, p);
    }
    return v;
}