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
    char type;
    uint64_t birthPos; // used in birth
    uint64_t moveDest; // used in move

    Atom *atom1; // used in death, move, exchange
    Atom *atom2; // used in exchange

    mutable GapsRng rng;

    // birth
    AtomicProposal(char t, uint64_t pos)
        : type(t), birthPos(pos), moveDest(0), atom1(NULL), atom2(NULL)
    {}
        
    // death
    AtomicProposal(char t, Atom *atom)
        : type(t), birthPos(0), moveDest(0), atom1(atom), atom2(NULL)
    {}

    // move
    AtomicProposal(char t, Atom *atom, uint64_t dest)
        : type(t), birthPos(0), moveDest(dest), atom1(atom), atom2(NULL)
    {}

    // exchange
    AtomicProposal(char t, Atom *a1, Atom *a2)
        : type(t), birthPos(0), moveDest(0), atom1(a1), atom2(a2)
    {}
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

    unsigned mNumBins; // number of matrix elements
    uint64_t mBinLength; // atomic length of one bin
    uint64_t mSecondaryDimLength; // atomic length of one row (col) for A (P)
    uint64_t mDomainLength; // length of entire atomic domain
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