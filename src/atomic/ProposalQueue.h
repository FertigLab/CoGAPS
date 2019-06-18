#ifndef __COGAPS_PROPOSAL_QUEUE_H__
#define __COGAPS_PROPOSAL_QUEUE_H__

#include "ConcurrentAtomicDomain.h"
#include "../data_structures/HashSets.h"
#include "../math/Random.h"
#include "../utils/Archive.h"

#include <cstddef>
#include <stdint.h>
#include <vector>

struct AtomicProposal
{
    mutable GapsRng rng; // used for consistency no matter number of threads 
 
    uint64_t pos; // used for move
    Atom *atom1; // used for birth/death/move/exchange
    Atom *atom2; // used for exchange

    uint32_t r1; // row of atom1
    uint32_t c1; // col of atom1
    uint32_t r2; // row of pos (move) or atom2 (exchange)
    uint32_t c2; // col of pos (move) or atom2 (exchange)

    char type; // birth (B), death (D), move (M), exchange (E)

    AtomicProposal(char t, GapsRandomState *randState);
};

class ProposalQueue
{
public:

    // initialize
    ProposalQueue(uint64_t nElements, uint64_t nPatterns, GapsRandomState *randState);
    void setAlpha(float alpha);
    void setLambda(float lambda);

    // modify/access queue
    void populate(ConcurrentAtomicDomain &domain, unsigned limit);
    void clear();
    unsigned size() const;
    AtomicProposal& operator[](int n);

    // update min/max atoms
    void acceptDeath();
    void rejectDeath();
    void acceptBirth();
    void rejectBirth();

    unsigned nProcessed() const;

    // serialization
    friend Archive& operator<<(Archive &ar, const ProposalQueue &queue);
    friend Archive& operator>>(Archive &ar, ProposalQueue &queue);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    std::vector<AtomicProposal> mQueue; // not really a queue for now
    
    FixedHashSetU32 mUsedMatrixIndices;
    SmallHashSetU64 mUsedAtoms;
    SmallPairedHashSetU64 mProposedMoves;

    GapsRandomState *mRandState;
    mutable GapsRng mRng;

    uint64_t mMinAtoms;
    uint64_t mMaxAtoms;
    uint64_t mBinLength; // length of single bin
    uint64_t mNumCols;

    double mAlpha;
    double mDomainLength; // length of entire atomic domain
    double mNumBins; // number of matrix elements

    float mLambda;
    float mU1;
    float mU2;
    
    unsigned mNumProcessed;

    bool mUseCachedRng;

    float deathProb(double nAtoms) const;

    bool makeProposal(ConcurrentAtomicDomain &domain);
    bool birth(ConcurrentAtomicDomain &domain);
    bool death(ConcurrentAtomicDomain &domain);
    bool move(ConcurrentAtomicDomain &domain);
    bool exchange(ConcurrentAtomicDomain &domain);
};

#endif // __COGAPS_PROPOSAL_QUEUE_H__