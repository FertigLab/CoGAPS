#ifndef __GAPS_PROPOSAL_QUEUE_H__
#define __GAPS_PROPOSAL_QUEUE_H__

#include "AtomicDomain.h"
#include "data_structures/HashSets.h"
#include "math/Random.h"
#include "utils/Archive.h"

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

    AtomicProposal(char t);
};

class ProposalQueue
{
public:

    // initialize
    ProposalQueue(unsigned nrow, unsigned ncol);
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

    // serialization
    friend Archive& operator<<(Archive &ar, ProposalQueue &queue);
    friend Archive& operator>>(Archive &ar, ProposalQueue &queue);

private:

    std::vector<AtomicProposal> mQueue; // not really a queue for now
    
    FixedHashSetU32 mUsedMatrixIndices;
    SmallHashSetU64 mUsedAtoms;

    mutable GapsRng mRng;

    uint64_t mMinAtoms;
    uint64_t mMaxAtoms;
    uint64_t mBinLength; // length of single bin
    uint64_t mNumCols;

    double mAlpha;
    double mDomainLength; // length of entire atomic domain
    double mNumBins; // number of matrix elements

    float mU1;
    float mU2;

    bool mUseCachedRng;

    float deathProb(double nAtoms) const;

    bool makeProposal(AtomicDomain &domain);
    bool birth(AtomicDomain &domain);
    bool death(AtomicDomain &domain);
    bool move(AtomicDomain &domain);
    bool exchange(AtomicDomain &domain);
};

#endif