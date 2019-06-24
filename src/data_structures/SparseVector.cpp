#include "SparseVector.h"
#include "Vector.h"
#include "../utils/Archive.h"
#include "../utils/GapsAssert.h"

// counts number of bits set below position
// internal expression clears all bits as high as pos or higher
// number of remaning bits is returned
static unsigned countLowerBits(uint64_t u, unsigned pos)
{
    return __builtin_popcountll(u & ((1ull << pos) - 1ull));
}

SparseVector::SparseVector(unsigned size)
    :
mSize(size),
mIndexBitFlags(size / 64 + 1, 0)
{}

SparseVector::SparseVector(const std::vector<float> &v)
    :
mSize(v.size()),
mIndexBitFlags(v.size() / 64 + 1, 0)
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (v[i] > 0.f)
        {
            mData.push_back(v[i]);
            mIndexBitFlags[i / 64] |= (1ull << (i % 64));
        }
    }
}

SparseVector::SparseVector(const Vector &v)
    :
mSize(v.size()),
mIndexBitFlags(v.size() / 64 + 1, 0)
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (v[i] > 0.f)
        {
            mData.push_back(v[i]);
            mIndexBitFlags[i / 64] |= (1ull << (i % 64));
        }
    }
}

unsigned SparseVector::size() const
{
    return mSize;
}

void SparseVector::insert(unsigned i, float v)
{
    GAPS_ASSERT(v > 0.f);
    GAPS_ASSERT(!(mIndexBitFlags[i / 64] & (1ull << (i % 64)))); // this data should not exist
    unsigned dataIndex = 0;
    for (unsigned j = 0; j < i / 64; ++j)
    {
        dataIndex += __builtin_popcountll(mIndexBitFlags[j]);
    }
    dataIndex += countLowerBits(mIndexBitFlags[i / 64], i % 64);
    mData.insert(mData.begin() + dataIndex, v);
    mIndexBitFlags[i / 64] |= (1ull << (i % 64));
}

Vector SparseVector::getDense() const
{
    Vector v(mSize);
    unsigned sparseNdx = 0;
    for (unsigned i = 0; i < mIndexBitFlags.size(); ++i)
    {
        uint64_t flags = mIndexBitFlags[i];
        while (flags != 0u)
        {
            unsigned ndx = __builtin_ffsll(flags) - 1;
            GAPS_ASSERT(sparseNdx < mData.size());
            v[64 * i + ndx] = mData[sparseNdx++];
            flags ^= 1ull << ndx;
        }
    }
    return v;
}

float SparseVector::at(unsigned n) const
{
    if ((mIndexBitFlags[n / 64] & (1ull << (n % 64))) == 0u)
    {
        return 0.f;
    }
    unsigned sparseNdx = 0;
    for (unsigned i = 0; i < n / 64; ++i)
    {
        sparseNdx += __builtin_popcountll(mIndexBitFlags[i]);
    }
    sparseNdx += countLowerBits(mIndexBitFlags[n / 64], n % 64);
    return mData[sparseNdx];
}

float SparseVector::getIthElement(unsigned n) const
{
    return mData[n];
}

unsigned SparseVector::nElements() const
{
    return mData.size();
}

Archive& operator<<(Archive &ar, const SparseVector &vec)
{
    ar << vec.mSize;
    for (unsigned i = 0; i < vec.mIndexBitFlags.size(); ++i)
    {
        ar << vec.mIndexBitFlags[i];
    }
    for (unsigned i = 0; i < vec.mData.size(); ++i)
    {
        ar << vec.mData[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, SparseVector &vec)
{
    unsigned sz = 0;
    ar >> sz;
    GAPS_ASSERT(sz == vec.mSize);

    for (unsigned i = 0; i < vec.mIndexBitFlags.size(); ++i)
    {
        ar >> vec.mIndexBitFlags[i];
    }
    for (unsigned i = 0; i < vec.mData.size(); ++i)
    {
        ar >> vec.mData[i];
    }
    return ar;
}
