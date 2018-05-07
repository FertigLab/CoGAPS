#ifndef __GAPS_PROPOSAL_QUEUE_H__
#define __GAPS_PROPOSAL_QUEUE_H__

#include "Archive.h"
#include "AtomicDomain.h"

#include <boost/unordered_set.hpp>
#include <stdint.h>
#include <cstddef>

struct AtomicProposal
{
    char type;
    uint64_t pos1;
    float mass1;
    uint64_t pos2;
    float mass2;

    AtomicProposal(char t, uint64_t p1, float m1=0.f, uint64_t p2=0, float m2=0.f)
        : type(t), pos1(p1), mass1(m1), pos2(p2), mass2(m2)
    {}
};

// generate single atomic proposal for now
class ProposalQueue
{
private:

    std::vector<AtomicProposal> mQueue; // not really a queue for now
    
    boost::unordered_set<uint64_t> mUsedIndices; // used rows/cols for A/P matrix
    std::set<uint64_t> mUsedPositions; // used positions in atomic domain

    uint64_t mMinAtoms;
    uint64_t mMaxAtoms;

    uint64_t mNumBins;
    uint64_t mDimensionSize; // rows of A, cols of P
    uint64_t mDomainSize;

    float mAlpha;

    float deathProb(uint64_t nAtoms) const;
    bool birth(AtomicDomain &domain);
    bool death(AtomicDomain &domain);
    bool move(AtomicDomain &domain);
    bool exchange(AtomicDomain &domain);

    bool makeProposal(AtomicDomain &domain);

public:

    ProposalQueue(unsigned nBins, float alpha)
        : mMinAtoms(0), mMaxAtoms(0), mNumBins(nBins), mAlpha(alpha)
    {}

    // set parameters
    void setNumBins(unsigned nBins);
    void setDomainSize(uint64_t size);
    void setAlpha(float alpha);
    void setDimensionSize(unsigned nIndices);

    // modify/access queue
    void populate(AtomicDomain &domain, unsigned limit);
    void clear(unsigned n);
    unsigned size() const;
    const AtomicProposal& operator[](int n) const;

    // update min/max atoms
    void acceptDeath();
    void rejectDeath();
    void rejectBirth();

    // serialization
    friend Archive& operator<<(Archive &ar, ProposalQueue &queue);
    friend Archive& operator>>(Archive &ar, ProposalQueue &queue);
};

#endif