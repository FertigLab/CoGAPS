#ifndef __COGAPS_SPARSE_ITERATOR_H__
#define __COGAPS_SPARSE_ITERATOR_H__

#include "HybridVector.h"
#include "SparseVector.h"

template <unsigned N, class Iter>
float get(const Iter &it);

// only allow this class to constructed with N=1,2,3
template <unsigned N>
class SparseIterator
{
private:
    SparseIterator() {}
};

template<>
class SparseIterator<1>
{
public:

    explicit SparseIterator(const SparseVector &v);

    bool atEnd() const;
    void next();
    unsigned getIndex() const;

private:

    friend float get<1>(const SparseIterator<1> &it);

    const SparseVector &mSparse;
    uint64_t mSparseFlags;
    unsigned mSparseIndex;
    unsigned mTotalIndices;
    unsigned mBigIndex;
    unsigned mSmallIndex;
    bool mAtEnd;
};

template<>
class SparseIterator<2>
{
public:

    SparseIterator(const SparseVector &v, const HybridVector &h);

    bool atEnd() const;
    void next();
    void calculateCommonFlags();
    void getFlags();
    unsigned getIndex() const;

private:

    template <class Iter>
    friend void gotoNextCommon(Iter &it);

    friend float get<1>(const SparseIterator<2> &it);
    friend float get<2>(const SparseIterator<2> &it);

    const SparseVector &mSparse;  
    const HybridVector &mHybrid;

    uint64_t mSparseFlags;
    uint64_t mHybridFlags;
    uint64_t mCommonFlags;

    unsigned mTotalIndices;
    unsigned mBigIndex;
    unsigned mSmallIndex;
    unsigned mSparseIndex;
    bool mAtEnd;
};

template<>
class SparseIterator<3>
{
public:

    SparseIterator(const SparseVector &v, const HybridVector &h1,
    const HybridVector &h2);

    bool atEnd() const;
    void next();
    void calculateCommonFlags();
    void getFlags();
    unsigned getIndex() const;

private:

    template <class Iter>
    friend void gotoNextCommon(Iter &it);

    friend float get<1>(const SparseIterator<3> &it);
    friend float get<2>(const SparseIterator<3> &it);
    friend float get<3>(const SparseIterator<3> &it);

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

#endif // __COGAPS_SPARSE_ITERATOR_H__