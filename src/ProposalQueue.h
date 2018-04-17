#ifndef __GAPS_PROPOSAL_QUEUE_H__
#define __GAPS_PROPOSAL_QUEUE_H__

#include "Archive.h"
#include "AtomicDomain.h"

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

    uint64_t mNumBins;
    uint64_t mDomainSize;
    float mAlpha;

    float deathProb(unsigned nAtoms) const;
    AtomicProposal birth(const AtomicDomain &domain);
    AtomicProposal death(const AtomicDomain &domain);
    AtomicProposal move(const AtomicDomain &domain);
    AtomicProposal exchange(const AtomicDomain &domain);

public:

    ProposalQueue(unsigned nBins, float alpha)
        : mNumBins(nBins), mAlpha(alpha)
    {}

    // set parameters
    void setNumBins(unsigned nBins);
    void setDomainSize(uint64_t size);
    void setAlpha(float alpha);

    // modify/access queue
    AtomicProposal makeProposal(const AtomicDomain &domain);

    // serialization
    friend Archive& operator<<(Archive &ar, ProposalQueue &queue);
    friend Archive& operator>>(Archive &ar, ProposalQueue &queue);
};

#endif