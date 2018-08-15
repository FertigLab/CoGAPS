#ifndef __GAPS_PROPOSAL_QUEUE_H__
#define __GAPS_PROPOSAL_QUEUE_H__

#include "Archive.h"
#include "AtomicDomain.h"
#include "data_structures/EfficientSets.h"
#include "math/Random.h"

#include <boost/unordered_set.hpp>

#include <cstddef>
#include <stdint.h>

struct AtomicProposal
{
    char type;
    uint64_t pos1;
    float mass1;
    uint64_t pos2;
    float mass2;
    mutable GapsRng rng;

    AtomicProposal(char t, uint64_t p1, float m1=0.f, uint64_t p2=0, float m2=0.f)
        : type(t), pos1(p1), mass1(m1), pos2(p2), mass2(m2)
    {}
};

class ProposalQueue
{
private:

    std::vector<AtomicProposal> mQueue; // not really a queue for now
    
    IntFixedHashSet mUsedIndices;
    IntDenseOrderedSet mUsedPositions;

    uint64_t mMinAtoms;
    uint64_t mMaxAtoms;

    uint64_t mNumBins;
    uint64_t mDimensionSize; // rows of A, cols of P
    uint64_t mDomainSize;

    float mAlpha;

    bool mUseCachedRng;
    float mU1;
    float mU2;

    mutable GapsRng mRng;

    float deathProb(uint64_t nAtoms) const;
    bool birth(AtomicDomain &domain);
    bool death(AtomicDomain &domain);
    bool move(AtomicDomain &domain);
    bool exchange(AtomicDomain &domain);

    bool makeProposal(AtomicDomain &domain);

public:

    ProposalQueue(unsigned nBins)
        : mMinAtoms(0), mMaxAtoms(0), mNumBins(nBins), mDimensionSize(0),
        mDomainSize(0), mAlpha(0.f), mUseCachedRng(false), mU1(0.f), mU2(0.f)
    {}

    // set parameters
    void setNumBins(unsigned nBins);
    void setDomainSize(uint64_t size);
    void setAlpha(float alpha);
    void setDimensionSize(uint64_t binSize, uint64_t dimLength);

    // modify/access queue
    void populate(AtomicDomain &domain, unsigned limit);
    void clear();
    unsigned size() const;
    const AtomicProposal& operator[](int n) const;

    // update min/max atoms
    void acceptDeath();
    void rejectDeath();
    void acceptBirth();
    void rejectBirth();

    // serialization
    friend Archive& operator<<(Archive &ar, ProposalQueue &queue);
    friend Archive& operator>>(Archive &ar, ProposalQueue &queue);
};

#endif