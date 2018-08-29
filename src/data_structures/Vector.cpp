#include "Vector.h"

#include "../utils/GapsAssert.h"

Vector::Vector(unsigned size)
    : mValues(aligned_vector(size, 0.f))
{}

Vector::Vector(const std::vector<float> &v) : mValues(v.size())
{
    unsigned sz = v.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        mValues[i] = v[i];
    }
}

unsigned Vector::size() const
{
    return mValues.size();
}

float* Vector::ptr()
{
    return &mValues[0];
}

const float* Vector::ptr() const
{
    return &mValues[0];
}

float& Vector::operator[](unsigned i)
{
    return mValues[i];
}

float Vector::operator[](unsigned i) const
{
    return mValues[i];
}

void Vector::operator+=(const Vector &vec)
{
    unsigned sz = size();
    for (unsigned i = 0; i < sz; ++i)
    {
        mValues[i] += vec[i];
    }
}

Vector Vector::operator-(Vector v) const
{
    unsigned sz = size();
    for (unsigned i = 0; i < sz; ++i)
    {
        v[i] = mValues[i] - v[i];
    }
    return v;
}

Vector Vector::operator*(float val) const
{
    Vector vec(*this);
    vec *= val;
    return vec;
}

Vector Vector::operator/(float val) const
{
    Vector vec(*this);
    vec /= val;
    return vec;
}

void Vector::operator*=(float val)
{
    unsigned sz = size();
    for (unsigned i = 0; i < sz; ++i)
    {
        mValues[i] *= val;
    }
}

void Vector::operator/=(float val)
{
    unsigned sz = size();
    for (unsigned i = 0; i < sz; ++i)
    {
        mValues[i] /= val;
    }
}

Archive& operator<<(Archive &ar, Vector &vec)
{
    unsigned sz = vec.size();
    ar << sz;
    for (unsigned i = 0; i < sz; ++i)
    {
        ar << vec[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, Vector &vec)
{
    unsigned sz = 0;
    ar >> sz;
    GAPS_ASSERT(sz == vec.size());

    for (unsigned i = 0; i < sz; ++i)
    {
        ar >> vec.mValues[i];
    }
    return ar;
}
