#ifndef __COGAPS_HASH_SETS_H__
#define __COGAPS_HASH_SETS_H__

#include <stdint.h>
#include <vector>

class FixedHashSetU32
{
public:
    explicit FixedHashSetU32(unsigned size);
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
private:
    std::vector<uint64_t> mSet;
};

struct PositionPair
{
    PositionPair(uint64_t inA, uint64_t inB) : a(inA), b(inB) {}
    uint64_t a;
    uint64_t b;
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
private:
    std::vector<PositionPair> mSet;
};

#endif // __COGAPS_HASH_SETS_H__
