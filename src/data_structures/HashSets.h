#ifndef __COGAPS_HASH_SETS_H__
#define __COGAPS_HASH_SETS_H__

#include "../utils/GlobalConfig.h"

#include <stdint.h>
#include <vector>

#ifdef __GAPS_OPENMP__ // defined in global config
#include <omp.h>
#endif

class FixedHashSetU32
{
public:

    explicit FixedHashSetU32(unsigned size);

    void insert(unsigned n);
    void clear();
    bool contains(unsigned n);
    bool isEmpty();

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

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

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    std::vector<uint64_t> mSet;
};

struct PositionPair
{
    uint64_t a;
    uint64_t b;

    PositionPair(uint64_t inA, uint64_t inB) : a(inA), b(inB) {}
};

class SmallPairedHashSetU64
{
public:

    SmallPairedHashSetU64();

    void insert(uint64_t a, uint64_t b);
    void clear();
    bool contains(uint64_t pos) const; // endpoint of pair
    bool overlap(uint64_t pos); // this position in between pair
    bool isEmpty();

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    std::vector<PositionPair> mSet;
};

#endif // __COGAPS_HASH_SETS_H__
