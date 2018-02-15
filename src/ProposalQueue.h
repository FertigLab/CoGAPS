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

// generate list of independent proposals
// does order really matter? i.e. could this be a stack?
// TODO need to cache failed proposal random numbers
// TODO use hash tables for atom conflicts
class ProposalQueue
{
private:

    std::vector<AtomicProposal> mQueue; // not really a queue for now

    std::vector<unsigned> mUsedIndices; // used rows/cols for A/P matrix
    std::vector<unsigned> mUsedPositions; // used positions is atomic domain

    uint64_t mMinAtoms;
    uint64_t mMaxAtoms;

    uint64_t mNumBins;
    uint64_t mNumIndices;

    float mAlpha;

    float deleteProb(unsigned nAtoms) const;
    bool makeProposal(const AtomicDomain &domain);
    bool birth(const AtomicDomain &domain);
    bool death(const AtomicDomain &domain);
    bool move(const AtomicDomain &domain);
    bool exchange(const AtomicDomain &domain);

public:

    // constructor
    ProposalQueue(float alpha, uint64_t nIndices, uint64_t nBins);

    // modify/access queue
    void populate(const AtomicDomain &domain, unsigned limit);
    void clear();
    unsigned size() const;
    const AtomicProposal& operator[](unsigned n) const;

    // update min/max atoms
    void acceptDeath();
    void rejectDeath();

    // serialization
    friend Archive& operator<<(Archive &ar, ProposalQueue &queue);
    friend Archive& operator>>(Archive &ar, ProposalQueue &queue);
};

#endif