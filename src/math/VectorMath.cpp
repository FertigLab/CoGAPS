#include "Math.h"
#include "VectorMath.h"
#include "../data_structures/SparseIterator.h"

float gaps::min(const Vector &v)
{
    float mn = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        mn = (v[i] < mn) ? v[i] : mn;
    }
    return mn;
}

float gaps::min(const HybridVector &v)
{
    float mn = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        mn = (v[i] < mn) ? v[i] : mn;
    }
    return mn;
}

float gaps::min(const SparseVector &v)
{
    float mn = 0.f;
    SparseIterator<1> it(v);
    while (!it.atEnd())
    {
        mn = (get<1>(it) < mn) ? get<1>(it) : mn;
        it.next();
    }
    return mn;
}

float gaps::max(const Vector &v)
{
    float mx = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        mx = (v[i] > mx) ? v[i] : mx;
    }
    return mx;
}

float gaps::max(const HybridVector &v)
{
    float mx = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        mx = (v[i] > mx) ? v[i] : mx;
    }
    return mx;
}

float gaps::max(const SparseVector &v)
{
    float mx = 0.f;
    SparseIterator<1> it(v);
    while (!it.atEnd())
    {
        mx = (get<1>(it) > mx) ? get<1>(it) : mx;
        it.next();
    }
    return mx;
}

unsigned gaps::whichMax(const Vector &v)
{
    unsigned ndx = 0;
    float mx = 0.f;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        ndx = (v[i] > mx) ? i : ndx;
        mx = (v[i] > mx) ? v[i] : mx;
    }
    return ndx;
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

float gaps::sum(const SparseVector &v)
{
    float sum = 0.f;
    SparseIterator<1> it(v);
    while (!it.atEnd())
    {
        sum += get<1>(it);
        it.next();
    }
    return sum;
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