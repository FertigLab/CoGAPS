#include "Vector.h"
#include "../math/SIMD.h"
#include "../utils/Archive.h"
#include "../utils/GapsAssert.h"

// Minimum effective SIMD width for padding purposes.
// When explicit SIMD is disabled (SIMD_INC == 1), the compiler may still
// auto-vectorize loops with 4-wide NEON on ARM or 4-wide SSE on x86.
// Use 8 to cover AVX2 as well.
#define SIMD_PAD_INC 8
#define SIMD_PAD(x) (SIMD_PAD_INC + SIMD_PAD_INC * ((x) / SIMD_PAD_INC))

Vector::Vector(unsigned sz)
    :
mData(SIMD_PAD(sz), 0.f),
mSize(sz)
{
    GAPS_ASSERT((mData.size() % gaps::simd::Index::increment()) == 0);
}

Vector::Vector(const std::vector<float> &v)
    :
mData(SIMD_PAD(v.size()), 0.f),
mSize(v.size())
{
    GAPS_ASSERT((mData.size() % gaps::simd::Index::increment()) == 0);
    for (unsigned i = 0; i < mSize; ++i)
    {
        mData[i] = v[i];
    }
}

void Vector::pad(float val)
{
    //here was the Msize as start value, error!
    for (unsigned i = 0; i < mData.size(); ++i)
    {
        mData[i] = val;
    }
}

// fill only the SIMD overflow positions beyond mSize with val
// to prevent 0/0 = NaN when SIMD reads past the real data
void Vector::padSIMD(float val)
{
    for (unsigned i = mSize; i < mData.size(); ++i)
    {
        mData[i] = val;
    }
}

float Vector::operator[](unsigned i) const
{
    GAPS_ASSERT(i < mSize);
    return mData[i];
}

float& Vector::operator[](unsigned i)
{
    GAPS_ASSERT(i < mSize);
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

void Vector::operator*=(float f)
{
    for (unsigned i = 0; i < mSize; ++i)
    {
        GAPS_ASSERT_MSG(mData[i] >= 0.f, mData[i]);
        mData[i] *= f;
    }
}

void Vector::operator/=(float f)
{
    for (unsigned i = 0; i < mSize; ++i)
    {
        GAPS_ASSERT_MSG(mData[i] >= 0.f, i << " , " << mSize << " : " << mData[i]);
        mData[i] /= f;
    }
}

Archive& operator<<(Archive &ar, const Vector &vec)
{
    ar << vec.mSize;
    for (unsigned i = 0; i < vec.mSize; ++i)
    {
        ar << vec.mData[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, Vector &vec)
{
    unsigned sz = 0;
    ar >> sz;
    GAPS_ASSERT(sz == vec.mSize);

    for (unsigned i = 0; i < vec.mSize; ++i)
    {
        ar >> vec.mData[i];
    }
    return ar;
}
