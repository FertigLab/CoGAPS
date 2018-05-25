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

Vector Vector::operator*(float val) const
{
    Vector vec(*this);
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] *= val;
    }
    return vec;
}

Vector Vector::operator/(float val) const
{
    Vector vec(*this);
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] /= val;
    }
    return vec;
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
