#ifndef __GAPS_PROPOSAL_QUEUE_H__
#define __GAPS_PROPOSAL_QUEUE_H__

struct AtomicProposal
{
    char type;
    uint64_t pos1;
    uint64_t pos2;
    
    AtomicProposal(char t, uint64_t p1, uint64_t p2)
        : type(t), pos1(p1), pos2(p2)
    {}
};

// generate list of independent proposals
class ProposalQueue
{
private:

public:

    ProposalQueue();

    void populate(unsigned limit); // add up to limit proposals to the queue
    void clear(unsigned n); // delete first n proposals
    const AtomicProposal& operator[](unsigned n); // get nth proposal
};

#endif