#ifndef __GAPS_PROPOSAL_QUEUE_H__
#define __GAPS_PROPOSAL_QUEUE_H__

#include "Archive.h"
#include "AtomicDomain.h"
#include "data_structures/EfficientSets.h"
#include "math/Random.h"

#include <cstddef>
#include <stdint.h>
#include <vector>

struct AtomicProposal
{
    uint64_t pos;
    Atom *atom1; // used in death, move, exchange
    Atom *atom2; // used in exchange

    mutable GapsRng rng;

    char type;

    AtomicProposal(char t, Atom *a); // birth/death
    AtomicProposal(char t, Atom *a, uint64_t p); // move
    AtomicProposal(char t, Atom *a1, Atom *a2); // exchange
};

class ProposalQueue
{
public:

    ProposalQueue(unsigned primaryDimSize, unsigned secondaryDimSize);
    void setAlpha(float alpha);

    // modify/access queue
    void populate(AtomicDomain &domain, unsigned limit);
    void clear();
    unsigned size() const;
    AtomicProposal& operator[](int n);

    // update min/max atoms
    void acceptDeath();
    void rejectDeath();
    void acceptBirth();
    void rejectBirth();

private:

    std::vector<AtomicProposal> mQueue; // not really a queue for now
    
    IntFixedHashSet mUsedIndices;
    IntDenseOrderedSet mUsedPositions;

    uint64_t mMinAtoms;
    uint64_t mMaxAtoms;

    double mNumBins; // number of matrix elements
    uint64_t mBinLength; // atomic length of one bin
    uint64_t mSecondaryDimLength; // atomic length of one row (col) for A (P)
    double mDomainLength; // length of entire atomic domain
    unsigned mSecondaryDimSize; // number of cols (rows) for A (P)

    float mAlpha;

    mutable GapsRng mRng;

    float mU1;
    float mU2;
    bool mUseCachedRng;

    unsigned primaryIndex(uint64_t pos) const;
    unsigned secondaryIndex(uint64_t pos) const;

    float deathProb(uint64_t nAtoms) const;
    bool birth(AtomicDomain &domain);
    bool death(AtomicDomain &domain);
    bool move(AtomicDomain &domain);
    bool exchange(AtomicDomain &domain);

    bool makeProposal(AtomicDomain &domain);

    // serialization
    friend Archive& operator<<(Archive &ar, ProposalQueue &queue);
    friend Archive& operator>>(Archive &ar, ProposalQueue &queue);
};

#endif