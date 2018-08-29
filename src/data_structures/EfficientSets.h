#ifndef __COGAPS_EFFICIENT_SETS_H__
#define __COGAPS_EFFICIENT_SETS_H__

#include <vector>
#include <stdint.h>

class IntFixedHashSet
{
public:

    IntFixedHashSet() : mCurrentKey(1) {}

    void setDimensionSize(unsigned size) {mSet.resize(size, 0);}
    void clear() {++mCurrentKey;}
    bool contains(unsigned n) {return mSet[n] == mCurrentKey;}
    void insert(unsigned n) {mSet[n] = mCurrentKey;}

private:

    std::vector<uint64_t> mSet;
    uint64_t mCurrentKey;
};

// TODO have sorted vector with at least some % of holes
// even distribute entries along it
// when shift happens, should be minimal
class IntDenseOrderedSet
{
public:

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

private:

    std::vector<uint64_t> mVec;
};

#endif // __COGAPS_EFFICIENT_SETS_H__