#ifndef __COGAPS_SPARSE_ITERATOR_H__
#define __COGAPS_SPARSE_ITERATOR_H__

#include "HybridVector.h"
#include "SparseVector.h"

// TODO make these nicer with templates - make sure no performance lost

class SparseIterator
{
public:

    SparseIterator(const SparseVector &v);

    bool atEnd() const;
    void next();
    float getValue() const;

private:

    const SparseVector &mSparse;
    unsigned mSparseIndex;
};

class SparseIteratorTwo
{
public:

    SparseIteratorTwo(const SparseVector &v1, const HybridVector &v2);

    bool atEnd() const;
    void next();

    float getValue_1() const;
    float getValue_2() const;
    unsigned getIndex() const;

private:

    const SparseVector &mSparse;
    const HybridVector &mHybrid;

    uint64_t mFlags_1;
    uint64_t mFlags_2;
    uint64_t mCommon;

    unsigned mTotalIndices;
    unsigned mBigIndex;
    unsigned mSmallIndex;
    unsigned mSparseIndex;

    bool mAtEnd;
};

class SparseIteratorThree
{
public:

    SparseIteratorThree(const SparseVector &v1, const HybridVector &v2,
    const HybridVector &v3);

    bool atEnd() const;
    void next();

    float getValue_1() const;
    float getValue_2() const;
    float getValue_3() const;
    unsigned getIndex() const;

private:

    const SparseVector &mSparse;
    const HybridVector &mHybrid_1;
    const HybridVector &mHybrid_2;
   
    uint64_t mFlags_1;
    uint64_t mFlags_2;
    uint64_t mFlags_3;
    uint64_t mCommon;

    unsigned mTotalIndices;
    unsigned mBigIndex;
    unsigned mSmallIndex;
    unsigned mSparseIndex;

    bool mAtEnd;
};

#endif // __COGAPS_SPARSE_ITERATOR_H__