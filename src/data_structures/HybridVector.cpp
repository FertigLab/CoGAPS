#include "HybridVector.h"
#include "../math/Math.h"
#include "../math/SIMD.h"

#define PAD_SIZE_FOR_SIMD(x) (gaps::simd::Index::increment() * (1 + ((x) - 1) / gaps::simd::Index::increment()))

HybridVector::HybridVector(unsigned size)
    :
mIndexBitFlags(size / 64 + 1, 0),
mData(PAD_SIZE_FOR_SIMD(size), 0.f),
mSize(size)
{}

HybridVector::HybridVector(const std::vector<float> &v)
    :
mIndexBitFlags(v.size() / 64 + 1, 0),
mData(PAD_SIZE_FOR_SIMD(v.size()), 0.f),
mSize(v.size())
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        mData[i] = v[i];
        if (v[i] > 0.f)
        {
            mIndexBitFlags[i / 64] ^= (1ull << (i % 64));
        }
    }
}

bool HybridVector::empty() const
{
    for (unsigned i = 0; i < mIndexBitFlags.size(); ++i)
    {
        if (mIndexBitFlags[i] != 0)
        {
            return false;
        }
    }
    return true;
}

unsigned HybridVector::size() const
{
    return mSize;
}

bool HybridVector::add(unsigned i, float v)
{
    if (mData[i] + v < gaps::epsilon)
    {
        mIndexBitFlags[i / 64] ^= (1ull << (i % 64));
        mData[i] = 0.f;
        return true;
    }
    else
    {
        mIndexBitFlags[i / 64] |= (1ull << (i % 64));
        mData[i] += v;
        return false;
    }
}

float HybridVector::operator[](unsigned i) const
{
    return mData[i];
}

const float* HybridVector::densePtr() const
{
    return &(mData[0]);
}

Archive& operator<<(Archive &ar, HybridVector &vec)
{
    ar << vec.mSize;
    for (unsigned i = 0; i < vec.mIndexBitFlags.size(); ++i)
    {
        ar << vec.mIndexBitFlags[i];
    }

    for (unsigned i = 0; i < vec.mSize; ++i)
    {
        ar << vec.mData[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, HybridVector &vec)
{
    unsigned sz = 0;
    ar >> sz;
    GAPS_ASSERT(sz == vec.size());

    for (unsigned i = 0; i < vec.mIndexBitFlags.size(); ++i)
    {
        ar >> vec.mIndexBitFlags[i];
    }

    for (unsigned i = 0; i < vec.mSize; ++i)
    {
        ar >> vec.mData[i];
    }
    return ar;
}


