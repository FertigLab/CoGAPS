#include "SparseIterator.h"
#include "../utils/GapsAssert.h"

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
    if (pos == 63)
    {
        return 0; // TODO understand this bug - should happen automatically
    }
    return u & ~((1ull << (pos + 1ull)) - 1ull);
}

template <class Iter>
void gotoNextCommon(Iter &it)
{
    // get the common indices in this chunk
    it.calculateCommonFlags();

    // if nothing common in this chunk, find a chunk that has common indices
    while (!it.mCommonFlags)
    {
        // first count how many sparse indices we are skipping
        it.mSparseIndex += __builtin_popcountll(it.mSparseFlags);
        
        // advance to next chunk
        if (++it.mBigIndex == it.mTotalIndices)
        {   
            it.mAtEnd = true;
            return;
        }

        // update the flags
        it.getFlags();
        it.calculateCommonFlags();
    }

    // must have at least one common value, this is our index
    it.mSmallIndex = __builtin_ffsll(it.mCommonFlags) - 1;

    // find the number of skipped entries in the sparse vector
    it.mSparseIndex += 1 + countLowerBits(it.mSparseFlags, it.mSmallIndex);

    // clear out all skipped indices and the current index from the bitflags
    it.mSparseFlags = clearLowerBits(it.mSparseFlags, it.mSmallIndex);
}

template<>
float get<1>(const SparseIterator<1> &it)
{
    return it.mSparse.getIthElement(it.mSparseIndex);
}

template<>
float get<1>(const SparseIterator<2> &it)
{
    return it.mSparse.getIthElement(it.mSparseIndex);
}

template<>
float get<1>(const SparseIterator<3> &it)
{
    return it.mSparse.getIthElement(it.mSparseIndex);
}

template<>
float get<2>(const SparseIterator<2> &it)
{
    GAPS_ASSERT(it.mHybrid[64 * it.mBigIndex + it.mSmallIndex] > 0.f);
    return it.mHybrid[64 * it.mBigIndex + it.mSmallIndex];
}

template<>
float get<2>(const SparseIterator<3> &it)
{
    return it.mHybrid_1[64 * it.mBigIndex + it.mSmallIndex];
}

template<>
float get<3>(const SparseIterator<3> &it)
{
    return it.mHybrid_2[64 * it.mBigIndex + it.mSmallIndex];
}

SparseIterator<1>::SparseIterator(const SparseVector &v)
:
mSparse(v),
mSparseFlags(v.mIndexBitFlags[0]),
mSparseIndex(0),
mTotalIndices(v.mIndexBitFlags.size()),
mBigIndex(0),
mSmallIndex(0),
mAtEnd(false)
{
    next();
    mSparseIndex -= 1; // next puts us at position 1, this resets to 0
}

bool SparseIterator<1>::atEnd() const
{
    return mAtEnd;
}

void SparseIterator<1>::next()
{
    ++mSparseIndex;
    while (mSparseFlags == 0u)
    {
        // advance to next chunk
        if (++mBigIndex == mTotalIndices)
        {   
            mAtEnd = true;
            return;
        }
        mSparseFlags = mSparse.mIndexBitFlags[mBigIndex];
    }
    mSmallIndex = __builtin_ffsll(mSparseFlags) - 1;
    mSparseFlags = clearLowerBits(mSparseFlags, mSmallIndex);
}

unsigned SparseIterator<1>::getIndex() const
{
    return 64 * mBigIndex + mSmallIndex;
}

SparseIterator<2>::SparseIterator(const SparseVector &v, const HybridVector &h)
    :
mSparse(v),
mHybrid(h),
mSparseFlags(v.mIndexBitFlags[0]),
mHybridFlags(h.mIndexBitFlags[0]),
mCommonFlags(v.mIndexBitFlags[0] & h.mIndexBitFlags[0]),
mTotalIndices(v.mIndexBitFlags.size()),
mBigIndex(0),
mSmallIndex(0),
mSparseIndex(0),
mAtEnd(false)
{
    GAPS_ASSERT(v.size() == h.size());

    next();
    mSparseIndex -= 1; // next puts us at position 1, this resets to 0
}

bool SparseIterator<2>::atEnd() const
{
    return mAtEnd;
}

void SparseIterator<2>::next()
{
    gotoNextCommon(*this);
}

void SparseIterator<2>::calculateCommonFlags()
{
    mCommonFlags = mSparseFlags & mHybridFlags;
}

void SparseIterator<2>::getFlags()
{
    mSparseFlags = mSparse.mIndexBitFlags[mBigIndex];
    mHybridFlags = mHybrid.mIndexBitFlags[mBigIndex];
}

unsigned SparseIterator<2>::getIndex() const
{
    return 64 * mBigIndex + mSmallIndex;
}

SparseIterator<3>::SparseIterator(const SparseVector &v,
const HybridVector &h1, const HybridVector &h2)
    :
mSparse(v),
mHybrid_1(h1),
mHybrid_2(h2),
mSparseFlags(v.mIndexBitFlags[0]),
mHybridFlags_1(h1.mIndexBitFlags[0]),
mHybridFlags_2(h2.mIndexBitFlags[0]),
mCommonFlags(v.mIndexBitFlags[0] & (h1.mIndexBitFlags[0] | h2.mIndexBitFlags[0])),
mTotalIndices(v.mIndexBitFlags.size()),
mBigIndex(0),
mSmallIndex(0),
mSparseIndex(0),
mAtEnd(false)
{
    GAPS_ASSERT(v.size() == h1.size());
    GAPS_ASSERT(h1.size() == h2.size());

    next();
    mSparseIndex -= 1;
}

bool SparseIterator<3>::atEnd() const
{
    return mAtEnd;
}

void SparseIterator<3>::next()
{
    gotoNextCommon(*this);
}

void SparseIterator<3>::calculateCommonFlags()
{
    mCommonFlags = mSparseFlags & (mHybridFlags_1 | mHybridFlags_2);
}

void SparseIterator<3>::getFlags()
{
    mSparseFlags = mSparse.mIndexBitFlags[mBigIndex];
    mHybridFlags_1 = mHybrid_1.mIndexBitFlags[mBigIndex];
    mHybridFlags_2 = mHybrid_2.mIndexBitFlags[mBigIndex];
}

unsigned SparseIterator<3>::getIndex() const
{
    return 64 * mBigIndex + mSmallIndex;
}
