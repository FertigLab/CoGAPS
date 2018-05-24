#ifndef __COGAPS_EFFICIENT_SETS_H__
#define __COGAPS_EFFICIENT_SETS_H__

#include <stdint.h>

class IntFixedHashSet
{
private:

    std::vector<uint64_t> mSet;
    uint64_t mCurrentKey;

public:

    IntFixedHashSet() : mCurrentKey(1) {}

    void setDimensionSize(uint64_t size) {mSet.resize(size, 0);}
    void clear() {++mCurrentKey;}
    bool count(uint64_t n) {return mSet[n] == mCurrentKey;}
    void insert(uint64_t n) {mSet[n] = mCurrentKey;}
};

// TODO have sorted vector with at least some % of holes
// even distribute entries along it
// when shift happens, should be minimal

class IntDenseOrderedSet
{
private:

    std::vector<uint64_t> mVec;

public:

    IntDenseOrderedSet() {}

    void insert(uint64_t p) {mVec.push_back(p);}
    void clear() {mVec.clear();}

    // inclusive of a and b, TODO improve performance
    bool isEmptyInterval(uint64_t a, uint64_t b)
    {
        for (unsigned i = 0; i < mVec.size(); ++i)
        {
            if (mVec[i] >= a && mVec[i] <= b)
            {
                return false;
            }
        }
        return true;
    }
};

#endif