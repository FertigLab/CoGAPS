#include "HashSets.h"

///////////////////////////// FixedHashSetU32 //////////////////////////////////

FixedHashSetU32::FixedHashSetU32(unsigned size)
    : mSet(std::vector<uint32_t>(size, 0)), mCurrentKey(1)
{}

void FixedHashSetU32::insert(unsigned n)
{
    mSet[n] = mCurrentKey;
}

void FixedHashSetU32::clear()
{
    ++mCurrentKey;
}

bool FixedHashSetU32::contains(unsigned n)
{
    return mSet[n] == mCurrentKey;
}

bool FixedHashSetU32::isEmpty()
{
    unsigned sz = mSet.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        if (mSet[i] == mCurrentKey)
        {
            return false;
        }
    }
    return true;
}

///////////////////////////// SmallHashSetU64 //////////////////////////////////

SmallHashSetU64::SmallHashSetU64() {}

void SmallHashSetU64::insert(uint64_t pos)
{
    mSet.push_back(pos);
}

void SmallHashSetU64::clear()
{
    mSet.clear();
}

bool SmallHashSetU64::contains(uint64_t pos)
{
    unsigned sz = mSet.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        if (mSet[i] == pos)
        {
            return true;
        }
    }
    return false;
}

bool SmallHashSetU64::isEmpty()
{
    return mSet.empty();
}

const std::vector<uint64_t>& SmallHashSetU64::vec()
{
    return mSet;
}
