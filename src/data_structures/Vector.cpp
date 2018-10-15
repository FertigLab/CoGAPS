#include "Vector.h"

Vector::Vector(unsigned size) : mData(size, 0.f) {}

Vector::Vector(const std::vector<float> &v) : mData(v.size(), 0.f)
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        mData[i] = v[i];
    }
}

float Vector::operator[](unsigned i) const
{
    return mData[i];
}

float& Vector::operator[](unsigned i)
{
    return mData[i];
}

float* Vector::ptr()
{
    return &(mData[0]);
}

const float* Vector::ptr() const
{
    return &(mData[0]);
}

unsigned Vector::size() const
{
    return mData.size();
}

void Vector::operator+=(const Vector &v)
{
    for (unsigned i = 0; i < mData.size(); ++i)
    {
        mData[i] += v[i];
    }
}

Archive& operator<<(Archive &ar, Vector &vec)
{
    for (unsigned i = 0; i < vec.mData.size(); ++i)
    {
        ar << vec.mData[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, Vector &vec)
{
    for (unsigned i = 0; i < vec.mData.size(); ++i)
    {
        ar >> vec.mData[i];
    }
    return ar;
}