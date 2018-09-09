#ifndef __COGAPS_HASH_SETS_H__
#define __COGAPS_HASH_SETS_H__

#include "../utils/GlobalConfig.h"

#include <stdint.h>
#include <vector>

#include <boost/unordered_set.hpp>

#ifdef __GAPS_OPENMP__
#include <omp.h>
#endif

class FixedHashSetU32
{
public:

    FixedHashSetU32(unsigned size);

    void insert(unsigned n);
    void clear();
    bool contains(unsigned n);
    bool isEmpty();

private:

    std::vector<uint32_t> mSet;
    uint64_t mCurrentKey;
};

class SmallHashSetU64
{
public:

    SmallHashSetU64();

    void insert(uint64_t pos);
    void clear();
    bool contains(uint64_t pos);
    bool isEmpty();
    const std::vector<uint64_t>& vec();

private:

    std::vector<uint64_t> mSet;
};

#endif // __COGAPS_HASH_SETS_H__
