#include "HybridVector.h"
#include "../math/Math.h"
#include "../math/SIMD.h"
#include "../utils/Archive.h"
#include "../utils/GapsAssert.h"

#define SIMD_PAD(x) (gaps::simd::Index::increment() + \
    gaps::simd::Index::increment() * ((x) / gaps::simd::Index::increment())) 

HybridVector::HybridVector(unsigned sz)
    :
mIndexBitFlags(sz / 64 + 1, 0),
mData(SIMD_PAD(sz), 0.f),
mSize(sz)
{
    GAPS_ASSERT(mData.size() % gaps::simd::Index::increment() == 0);
}

HybridVector::HybridVector(const std::vector<float> &v)
    :
mIndexBitFlags(v.size() / 64 + 1, 0),
mData(SIMD_PAD(v.size()), 0.f),
mSize(v.size())
{
    GAPS_ASSERT(mData.size() % gaps::simd::Index::increment() == 0);

    for (unsigned i = 0; i < mSize; ++i)
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

// can be called from multiple concurrent OpenMP threads
bool HybridVector::add(unsigned i, float v)
{
    GAPS_ASSERT(i < mSize);
    if (mData[i] + v < gaps::epsilon)
    {
        #pragma omp atomic
        mIndexBitFlags[i / 64] &= ~(1ull << (i % 64));
        mData[i] = 0.f;
        return true;
    }
    #pragma omp atomic
    mIndexBitFlags[i / 64] |= (1ull << (i % 64));
    mData[i] += v;
    return false;
}

// can be called from multiple concurrent OpenMP threads
bool HybridVector::set(unsigned i, float v)
{
    GAPS_ASSERT(i < mSize);
    if (v < gaps::epsilon)
    {
        #pragma omp atomic
        mIndexBitFlags[i / 64] &= ~(1ull << (i % 64));
        mData[i] = 0.f;
        return true;
    }
    #pragma omp atomic
    mIndexBitFlags[i / 64] |= (1ull << (i % 64));
    mData[i] = v;
    return false;
}

float HybridVector::operator[](unsigned i) const
{
    GAPS_ASSERT(i < mSize);
    GAPS_ASSERT((mData[i] == 0.f)
        ? !(mIndexBitFlags[i / 64] & (1ull << (i % 64)))
        : (mIndexBitFlags[i / 64] & (1ull << (i % 64)))
    );
    return mData[i];
}

const float* HybridVector::ptr() const
{
    return &(mData[0]);
}

Archive& operator<<(Archive &ar, const HybridVector &vec)
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


