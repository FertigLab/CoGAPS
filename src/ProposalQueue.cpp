#include "GapsAssert.h"
#include "ProposalQueue.h"
#include "Random.h"

void ProposalQueue::setNumBins(unsigned nBins)
{
    mNumBins = nBins;
}

void ProposalQueue::setDomainSize(uint64_t size)
{
    mDomainSize = size;
}

void ProposalQueue::setAlpha(float alpha)
{
    mAlpha = alpha;
}

float ProposalQueue::deathProb(uint64_t nAtoms) const
{
    //double size = static_cast<double>(mDomainSize);
    //double term1 = (size - static_cast<double>(nAtoms)) / size;
    //double term2 = mAlpha * static_cast<double>(mNumBins) * term1;
    //return static_cast<double>(nAtoms) / (static_cast<double>(nAtoms) + term2);
    float dMax = (float)mDomainSize;
    float dNum = (float)nAtoms;
    float maxTerm = (dMax - dNum) / dMax;

    return dNum / (dNum + mAlpha * (float)mNumBins * maxTerm);
}

AtomicProposal ProposalQueue::makeProposal(const AtomicDomain &domain)
{
    // always birth when no atoms exist
    if (domain.size() == 0)
    {
        return birth(domain);
    }

    float bdProb = domain.size() < 2 ? 0.6667f : 0.5f;
    float u = gaps::random::uniform();
    if (u <= bdProb)
    {
        return gaps::random::uniform() < deathProb(domain.size()) ?
            death(domain) : birth(domain);
    }
    else if (u < 0.75f || domain.size() < 2)
    {
        return move(domain);
    }
    else
    {
        return exchange(domain);
    }
}

AtomicProposal ProposalQueue::birth(const AtomicDomain &domain)
{
    uint64_t pos = domain.randomFreePosition();
    return AtomicProposal('B', pos);
}

AtomicProposal ProposalQueue::death(const AtomicDomain &domain)
{
    Atom a = domain.randomAtom();
    GAPS_ASSERT(domain.test(a.pos));
    return AtomicProposal('D', a.pos, a.mass);
}

AtomicProposal ProposalQueue::move(const AtomicDomain &domain)
{
    Atom a = domain.randomAtom();
    uint64_t lbound = domain.hasLeft(a) ? domain.left(a).pos : 0;
    uint64_t rbound = domain.hasRight(a) ? domain.right(a).pos : mDomainSize;
    uint64_t newLocation = gaps::random::uniform64(lbound, rbound - 1);
    GAPS_ASSERT(domain.test(a.pos));
    return AtomicProposal('M', a.pos, a.mass, newLocation);
}

AtomicProposal ProposalQueue::exchange(const AtomicDomain &domain)
{
    Atom a1 = domain.randomAtom();
    Atom a2 = domain.hasRight(a1) ? domain.right(a1) : domain.front();
    GAPS_ASSERT(domain.test(a1.pos));
    GAPS_ASSERT(domain.test(a2.pos));
    GAPS_ASSERT(a1.pos != a2.pos);
    return AtomicProposal('E', a1.pos, a1.mass, a2.pos, a2.mass);
}