#include "Math.h"
#include "VectorMath.h"
#include "../data_structures/SparseIterator.h"

float gaps::min(const Vector &v)
{
    if (v.size() == 0) return 0.f; // empty vector
    float mn = v[0];
    for (unsigned i = 0; i < v.size(); ++i)
    {
        // [AI-generated] Keep the smaller value seen so far.
        mn = (v[i] < mn) ? v[i] : mn;
    }
    return mn;
}

float gaps::min(const HybridVector &v)
{
    if (v.size() == 0) return 0.f; // empty vector
    float mn = v[0];
    for (unsigned i = 0; i < v.size(); ++i)
    {
        // [AI-generated] Keep the smaller value seen so far.
        mn = (v[i] < mn) ? v[i] : mn;
    }
    return mn;
}

float gaps::min(const SparseVector &v)
{
    SparseIterator<1> it(v);
    if (it.atEnd()) return 0.f; // empty (all-zero) sparse vector
    float mn = get<1>(it);
    while (!it.atEnd())
    {
        // [AI-generated] Keep the smaller nonzero sparse value seen so far.
        mn = (get<1>(it) < mn) ? get<1>(it) : mn;
        it.next();
    }
    return mn;
}

float gaps::max(const Vector &v)
{
    if (v.size() == 0) return 0.f; // empty vector
    float mx = v[0];
    for (unsigned i = 1; i < v.size(); ++i)
    {
        // [AI-generated] Keep the larger value seen so far.
        mx = (v[i] > mx) ? v[i] : mx;
    }
    return mx;
}

float gaps::max(const HybridVector &v)
{
    if (v.size() == 0) return 0.f; // empty vector
    float mx = v[0];
    for (unsigned i = 1; i < v.size(); ++i)
    {
        // [AI-generated] Keep the larger value seen so far.
        mx = (v[i] > mx) ? v[i] : mx;
    }
    return mx;
}

float gaps::max(const SparseVector &v)
{
    SparseIterator<1> it(v);
    if (it.atEnd()) return 0.f; // empty (all-zero) sparse vector
    float mx = get<1>(it);
    while (!it.atEnd())
    {
        // [AI-generated] Keep the larger nonzero sparse value seen so far.
        mx = (get<1>(it) > mx) ? get<1>(it) : mx;
        it.next();
    }
    return mx;
}

unsigned gaps::whichMax(const Vector &v)
{
    if (v.size() == 0) return 0; // empty vector
    unsigned ndx = 0;
    float mx = v[0];
    for (unsigned i = 0; i < v.size(); ++i)
    {
        // [AI-generated] Track the index and value of the largest entry seen so far.
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

Vector gaps::elementSq(const Vector & v)
{
    Vector res(v.size());
    for (unsigned i = 0; i < v.size(); ++i)
    {
        res[i] = v[i] * v[i];
    }
    return res;
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

Vector gaps::pmax(const Vector & v, float f, float min_thr)
{
    Vector res(v.size());
    for (unsigned i = 0; i < v.size(); ++i)
    {
        res[i] = std::max(v[i] * f, min_thr);
    }
    res.padSIMD(min_thr);
    return res;
}

Vector gaps::pmax(const Vector & v, float f) {
    return (gaps::pmax(v,f,f));
}
