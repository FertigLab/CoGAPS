#include "Vector.h"
#include "../math/SIMD.h"

#define PAD_SIZE_FOR_SIMD(x) (gaps::simd::Index::increment() * (1 + ((x) - 1) / gaps::simd::Index::increment()))

Vector::Vector(unsigned size)
    :
mData(PAD_SIZE_FOR_SIMD(size), 0.f),
mSize(size)
{}

Vector::Vector(const std::vector<float> &v)
    :
mData(PAD_SIZE_FOR_SIMD(v.size()), 0.f),
mSize(v.size())
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
    return mSize;
}

void Vector::operator+=(const Vector &v)
{
    for (unsigned i = 0; i < mSize; ++i)
    {
        mData[i] += v[i];
    }
}

Archive& operator<<(Archive &ar, Vector &vec)
{
    for (unsigned i = 0; i < vec.mSize; ++i)
    {
        ar << vec.mData[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, Vector &vec)
{
    for (unsigned i = 0; i < vec.mSize; ++i)
    {
        ar >> vec.mData[i];
    }
    return ar;
}