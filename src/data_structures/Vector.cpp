#include "Vector.h"

Vector::Vector(const std::vector<float> &v) : mValues(v.size())
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        mValues[i] = v[i];
    }
}

void Vector::concat(const Vector& vec)
{
    mValues.insert(mValues.end(), vec.mValues.begin(), vec.mValues.end());
}

void Vector::operator+=(const Vector &vec)
{
    for (unsigned i = 0; i < size(); ++i)
    {
        mValues[i] += vec[i];
    }
}

Vector Vector::operator-(Vector v) const
{
    for (unsigned i = 0; i < size(); ++i)
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
    for (unsigned i = 0; i < mValues.size(); ++i)
    {
        mValues[i] *= val;
    }
}

void Vector::operator/=(float val)
{
    for (unsigned i = 0; i < size(); ++i)
    {
        mValues[i] /= val;
    }
}

Archive& operator<<(Archive &ar, Vector &vec)
{
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        ar << vec[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, Vector &vec)
{
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        ar >> vec.mValues[i];
    }
    return ar;
}
