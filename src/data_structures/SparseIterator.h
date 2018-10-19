#ifndef __COGAPS_SPARSE_ITERATOR_H__
#define __COGAPS_SPARSE_ITERATOR_H__

#include "HybridVector.h"
#include "SparseVector.h"

template <unsigned N, class Iter>
float get(const Iter &it);

// only allow this class to constructed with N=1,2,3
template <unsigned N>
class TemplatedSparseIterator
{
private:
    TemplatedSparseIterator() {}
};

template<>
class TemplatedSparseIterator<1>
{
public:

    TemplatedSparseIterator(const SparseVector &v);

    bool atEnd() const;
    void next();

private:

    friend float get<1>(const TemplatedSparseIterator<1> &it);

    const SparseVector &mSparse;
    unsigned mSparseIndex;
};

template<>
class TemplatedSparseIterator<2>
{
public:

    TemplatedSparseIterator(const SparseVector &v, const HybridVector &h);

    bool atEnd() const;
    void next();
    void calculateCommonFlags();
    void getFlags();

private:

    template <class Iter>
    friend void gotoNextCommon(Iter &it);

    friend float get<1>(const TemplatedSparseIterator<2> &it);
    friend float get<2>(const TemplatedSparseIterator<2> &it);

    const SparseVector &mSparse;  
    const HybridVector &mHybrid_1;

    uint64_t mSparseFlags;
    uint64_t mHybridFlags_1;
    uint64_t mCommonFlags;

    unsigned mTotalIndices;
    unsigned mBigIndex;
    unsigned mSmallIndex;
    unsigned mSparseIndex;
    bool mAtEnd;
};

template<>
class TemplatedSparseIterator<3>
{
public:

    TemplatedSparseIterator(const SparseVector &v, const HybridVector &h1,
    const HybridVector &h2);

    bool atEnd() const;
    void next();
    void calculateCommonFlags();
    void getFlags();

private:

    template <class Iter>
    friend void gotoNextCommon(Iter &it);

    friend float get<1>(const TemplatedSparseIterator<3> &it);
    friend float get<2>(const TemplatedSparseIterator<3> &it);
    friend float get<3>(const TemplatedSparseIterator<3> &it);

    const SparseVector &mSparse;  
    const HybridVector &mHybrid_1;
    const HybridVector &mHybrid_2;

    uint64_t mSparseFlags;
    uint64_t mHybridFlags_1;
    uint64_t mHybridFlags_2;
    uint64_t mCommonFlags;

    unsigned mTotalIndices;
    unsigned mBigIndex;
    unsigned mSmallIndex;
    unsigned mSparseIndex;
    bool mAtEnd;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

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