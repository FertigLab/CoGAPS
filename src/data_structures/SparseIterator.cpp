#include "SparseIterator.h"

// counts number of bits set below position
// internal expression clears all bits as high as pos or higher
// number of remaning bits is returned
static unsigned countLowerBits(uint64_t u, unsigned pos)
{
    return __builtin_popcountll(u & ((1ull << pos) - 1ull));
}

// clears all bits as low as pos or lower
static uint64_t clearLowerBits(uint64_t u, unsigned pos)
{
    return u & ~((1ull << (pos + 1ull)) - 1ull);    
}

SparseIterator::SparseIterator(const SparseVector &v)
    :
mSparse(v),
mSparseIndex(0)
{}

bool SparseIterator::atEnd() const
{
    return mSparseIndex == mSparse.mData.size();
}

void SparseIterator::next()
{
    ++mSparseIndex;
}

float SparseIterator::getValue() const
{
    return mSparse.mData[mSparseIndex];
}

SparseIteratorTwo::SparseIteratorTwo(const SparseVector &v1, const HybridVector &v2)
    :
mSparse(v1),
mHybrid(v2),
mFlags_1(mSparse.mIndexBitFlags[0]),
mFlags_2(mHybrid.mIndexBitFlags[0]),
mCommon(mFlags_1 & mFlags_2),
mTotalIndices(mSparse.mIndexBitFlags.size()),
mBigIndex(0),
mSmallIndex(0),
mSparseIndex(0),
mAtEnd(false)
{   
    GAPS_ASSERT(v1.size() == v2.size());
    next();
    mSparseIndex -= 1; // this gets advanced one too far
}

bool SparseIteratorTwo::atEnd() const
{
    return mAtEnd;
}

void SparseIteratorTwo::next()
{
    while (!mCommon)
    {
        // no common values in this chunk, go to next one
        ++mBigIndex;
        if (mBigIndex == mTotalIndices)
        {   
            mAtEnd = true;
            return;
        }

        // check if we have any common values here
        mFlags_1 = mSparse.mIndexBitFlags[mBigIndex];
        mFlags_2 = mHybrid.mIndexBitFlags[mBigIndex];
        mCommon = mFlags_1 & mFlags_2;
    }

    // must have at least one common value, this is our index
    mSmallIndex = __builtin_ffsll(mCommon) - 1;

    // find the number of skipped entries in the sparse vector
    mSparseIndex += 1 + countLowerBits(mFlags_1, mSmallIndex);

    // clear out all skipped indices and the current index from the bitflags
    // this is needed so that countLowerBits is accurate in the next iteration
    mFlags_1 = clearLowerBits(mFlags_1, mSmallIndex);

    // clear out this bit from common
    mCommon ^= (1ull << mSmallIndex);
}

float SparseIteratorTwo::getValue_1() const
{
    return mSparse.mData[mSparseIndex];
}

float SparseIteratorTwo::getValue_2() const
{
    return mHybrid[64 * mBigIndex + mSmallIndex];
}

unsigned SparseIteratorTwo::getIndex() const
{
    return 64 * mBigIndex + mSmallIndex;
}

SparseIteratorThree::SparseIteratorThree(const SparseVector &v1,
const HybridVector &v2, const HybridVector &v3)
    :
mSparse(v1),
mHybrid_1(v2),
mHybrid_2(v3),
mFlags_1(mSparse.mIndexBitFlags[0]),
mFlags_2(mHybrid_1.mIndexBitFlags[0]),
mFlags_3(mHybrid_2.mIndexBitFlags[0]),
mCommon(mFlags_1 & mFlags_2 & mFlags_3),
mTotalIndices(mSparse.mIndexBitFlags.size()),
mBigIndex(0),
mSmallIndex(0),
mSparseIndex(0),
mAtEnd(false)
{   
    GAPS_ASSERT(v1.size() == v2.size());
    GAPS_ASSERT(v2.size() == v3.size());
    next();
    mSparseIndex -= 1; // this gets advanced one too far
}

bool SparseIteratorThree::atEnd() const
{
    return mAtEnd;
}

void SparseIteratorThree::next()
{
    while (!mCommon)
    {
        // no common values in this chunk, go to next one
        ++mBigIndex;
        if (mBigIndex == mTotalIndices)
        {   
            mAtEnd = true;
            return;
        }

        // check if we have any common values here
        mFlags_1 = mSparse.mIndexBitFlags[mBigIndex];
        mFlags_2 = mHybrid_1.mIndexBitFlags[mBigIndex];
        mFlags_3 = mHybrid_2.mIndexBitFlags[mBigIndex];
        mCommon = mFlags_1 & mFlags_2 & mFlags_3;
    }

    // must have at least one common value, this is our index
    mSmallIndex = __builtin_ffsll(mCommon) - 1;

    // find the number of skipped entries in the sparse vector
    mSparseIndex += 1 + countLowerBits(mFlags_1, mSmallIndex);

    // clear out all skipped indices and the current index from the bitflags
    // this is needed so that countLowerBits is accurate in the next iteration
    mFlags_1 = clearLowerBits(mFlags_1, mSmallIndex);

    // clear out this bit from common
    mCommon ^= (1ull << mSmallIndex);
}

float SparseIteratorThree::getValue_1() const
{
    return mSparse.mData[mSparseIndex];
}

float SparseIteratorThree::getValue_2() const
{
    return mHybrid_1[64 * mBigIndex + mSmallIndex];
}

float SparseIteratorThree::getValue_3() const
{
    return mHybrid_2[64 * mBigIndex + mSmallIndex];
}

unsigned SparseIteratorThree::getIndex() const
{
    return 64 * mBigIndex + mSmallIndex;
}
