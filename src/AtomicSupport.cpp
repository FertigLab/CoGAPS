#include "AtomicSupport.h"

static const float EPSILON = 1.e-10;

AtomicSupport::AtomicSupport(char label, uint64_t nrow, uint64_t ncol,
float alpha, float lambda)
    :
mLabel(label), mNumAtoms(0),
mMaxNumAtoms(std::numeric_limits<uint64_t>::max()),
mNumRows(nrow), mNumCols(ncol), mNumBins(nrow * ncol),
mBinSize(std::numeric_limits<uint64_t>::max() / (nrow * ncol)),
mAlpha(alpha), mLambda(lambda)
{}

// O(logN)
void AtomicSupport::addAtom(Atom atom)
{
    mAtomPositions.insert(std::pair<uint64_t, uint64_t>(atom.pos, mNumAtoms));
    mAtoms.push_back(atom);
    mNumAtoms++;
}

// O(logN)
void AtomicSupport::removeAtom(uint64_t pos)
{
    uint64_t index = mAtomPositions.at(pos);
    mAtomPositions.erase(pos);

    // update key of object about to be moved (last one doesn't need to move)
    if (index < mNumAtoms - 1)
    {
        mAtomPositions.erase(mAtoms.back().pos);
        mAtomPositions.insert(std::pair<uint64_t, uint64_t>(mAtoms.back().pos,
            index));
    }

    // delete atom from vector in O(1)
    mAtoms[index] = mAtoms.back();
    mAtoms.pop_back();
    mNumAtoms--;
}

// O(logN)
AtomNeighbors AtomicSupport::getNeighbors(uint64_t pos) const
{
    // get an iterator to this atom and it's neighbors
    std::map<uint64_t, uint64_t>::const_iterator it, left, right;
    it = mAtomPositions.find(pos);
    left = it; right = it;
    if (it != mAtomPositions.begin())
    {
        --left;
    }
    if (++it != mAtomPositions.end())
    {
        ++right;
    }
    return AtomNeighbors(left->first, mAtoms[left->second].mass,
        right->first, mAtoms[right->second].mass);
}

// O(1)
uint64_t AtomicSupport::getRow(uint64_t pos) const
{
    uint64_t binNum = std::min(pos / mBinSize, mNumBins - 1);
    return mLabel == 'A' ? binNum / mNumCols : binNum % mNumRows;
}

// O(1)
uint64_t AtomicSupport::getCol(uint64_t pos) const
{
    uint64_t binNum = std::min(pos / mBinSize, mNumBins - 1);
    return mLabel == 'A' ? binNum % mNumCols : binNum / mNumRows;
}

// O(logN)
uint64_t AtomicSupport::randomFreePosition() const
{
    uint64_t pos = 0;
    do
    {
        pos = gaps::random::uniform64();
    } while (mAtomPositions.count(pos) > 0);
    return pos;
}

// O(1)
uint64_t AtomicSupport::randomAtomPosition() const
{
    uint64_t num = gaps::random::uniform64(0, mNumAtoms - 1);
    return mAtoms[num].pos;
}

// O(logN)
AtomicProposal AtomicSupport::proposeDeath() const
{
    uint64_t location = randomAtomPosition();
    float mass = mAtoms[mAtomPositions.at(location)].mass;
    return AtomicProposal(mLabel, 'D', location, -mass);
}

// O(logN)
AtomicProposal AtomicSupport::proposeBirth() const
{
    uint64_t location = randomFreePosition();
    float mass = gaps::random::exponential(mLambda);
    return AtomicProposal(mLabel, 'B', location, mass);
}

// move atom between adjacent atoms
// O(logN)
AtomicProposal AtomicSupport::proposeMove() const
{
    uint64_t pos = randomAtomPosition();
    AtomNeighbors nbs = getNeighbors(pos);
    float mass = at(pos);

    uint64_t lbound = nbs.left.pos != pos ? nbs.left.pos : 0;
    uint64_t rbound = nbs.right.pos != pos ? nbs.right.pos : mMaxNumAtoms - 1;

    uint64_t newLocation = gaps::random::uniform64(lbound, rbound);
    return AtomicProposal(mLabel, 'M', pos, -mass, newLocation, mass);
}

// exchange with adjacent atom to the right
// O(logN)
AtomicProposal AtomicSupport::proposeExchange() const
{
    uint64_t pos1 = randomAtomPosition();
    AtomNeighbors nbs = getNeighbors(pos1);
    float mass1 = at(pos1);

    // find right neighbor, wrap around to beginning if this atom is end
    uint64_t pos2 = nbs.right.pos != pos1 ? nbs.right.pos : mAtomPositions.begin()->first;
    float mass2 = at(pos2);

    // calculate new mass
    float pupper = gaps::random::p_gamma(mass1 + mass2, 2.f, 1.f / mLambda);
    float newMass = gaps::random::q_gamma(gaps::random::uniform(0.f, pupper),
        2.f, 1.f / mLambda);

    // calculate mass changes
    float delta1 = mass1 > mass2 ? newMass - mass1 : mass2 - newMass;
    float delta2 = mass2 > mass1 ? newMass - mass2 : mass1 - newMass;

    // set first position to be postive change
    float delta1_cached = delta1;
    uint64_t pos1_cached = pos1;
    pos1 = delta1_cached > 0 ? pos1 : pos2;
    pos2 = delta1_cached > 0 ? pos2 : pos1_cached;
    delta1 = delta1_cached > 0 ? delta1 : delta2;
    delta2 = delta1_cached > 0 ? delta2 : delta1_cached;

    // preserve total mass
    return AtomicProposal(mLabel, 'E', pos1, delta1, pos2, delta2);
}

// O(logN)
float AtomicSupport::updateAtomMass(uint64_t pos, float delta)
{
    if (mAtomPositions.count(pos)) // update atom if it exists
    {
        mAtoms[mAtomPositions.at(pos)].mass += delta;
    }
    else // create a new atom if it doesn't exist
    {
        addAtom(Atom(pos, delta));
    }

    if (mAtoms[mAtomPositions.at(pos)].mass < EPSILON)
    {
        delta -= at(pos);
        removeAtom(pos);
    }
    return delta;
}

MatrixChange AtomicSupport::acceptProposal(const AtomicProposal &prop,
MatrixChange &ch)
{
    ch.delta1 = updateAtomMass(prop.pos1, prop.delta1);
    ch.delta2 = (prop.nChanges > 1) ? updateAtomMass(prop.pos2, prop.delta2)
        : ch.delta2;
    return ch;
}

float AtomicSupport::deleteProb(unsigned nAtoms) const
{
    float pDelete = 0.5; // default uniform prior for mAlpha == 0
    float fnAtoms = (float)nAtoms;
    if (mAlpha > 0) // poisson prior
    {
        float maxTerm = ((float)mMaxNumAtoms - fnAtoms) / (float)mMaxNumAtoms;
        return fnAtoms / (fnAtoms + mAlpha * (float)mNumBins * maxTerm);
    }
    else if (mAlpha < 0) // geometric prior
    {
        float c = -mAlpha * mNumBins / (-mAlpha * mNumBins + 1.0);
        return fnAtoms / ((fnAtoms + 1.f) * c + fnAtoms);
    }    
}

bool AtomicSupport::addProposal(unsigned &minAtoms, unsigned &maxAtoms) const
{
    if (minAtoms == 0 && maxAtoms > 0 || minAtoms < 2 && maxAtoms >= 2)
    {
        return false;
    }

    if (maxAtoms == 0)
    {
        mProposalQueue.push_back(proposeBirth());
        minAtoms++;
        maxAtoms++;
        return true;
    }

    float bdProb = maxAtoms < 2 ? 0.6667f : 0.5f;
    float u1 = gaps::random::uniform();
    float u2 = gaps::random::uniform();   
    float lowerBound = deleteProb(minAtoms);
    float upperBound = deleteProb(maxAtoms);
    if (u1 < bdProb && u2 > upperBound)
    {
        mProposalQueue.push_back(proposeBirth());
        minAtoms++;
        maxAtoms++;
    }
    else if (u1 < bdProb && u2 < lowerBound)
    {
        mProposalQueue.push_back(proposeDeath());
        minAtoms--;
    }
    else if (u1 >= bdProb && (maxAtoms < 2 || u1 < 0.75f))
    {
        mProposalQueue.push_back(proposeMove());
    }
    else if (u1 >= bdProb)
    {
        mProposalQueue.push_back(proposeExchange());
        minAtoms--;
    }
    else
    {
        return false;
    }
    return true;
}

void AtomicSupport::populateQueue(unsigned limit)
{
    unsigned minAtoms = mNumAtoms;
    unsigned maxAtoms = mNumAtoms;
    unsigned nIter = 0;

    while (nIter++ < limit && addProposal(minAtoms, maxAtoms));
}

bool AtomicSupport::isQueueEmpty() const
{
    return mProposalQueue.empty();
}

AtomicProposal AtomicSupport::popQueue()
{
    AtomicProposal prop = mProposalQueue.back();
    mProposalQueue.pop_back();
}

Archive& operator<<(Archive &ar, AtomicSupport &domain)
{
    ar << domain.mLabel << domain.mNumAtoms << domain.mMaxNumAtoms
        << domain.mNumRows << domain.mNumCols
        << domain.mNumBins << domain.mBinSize << domain.mAlpha
        << domain.mLambda;

    for (unsigned i = 0; i < domain.mAtoms.size(); ++i)
    {
        ar << domain.mAtoms[i].pos << domain.mAtoms[i].mass;
    }
    return ar;
}

Archive& operator>>(Archive &ar, AtomicSupport &domain)
{
    uint64_t nAtoms = 0;
    ar >> domain.mLabel >> nAtoms >> domain.mMaxNumAtoms
        >> domain.mNumRows >> domain.mNumCols >> domain.mNumBins
        >> domain.mBinSize >> domain.mAlpha >> domain.mLambda;

    uint64_t pos = 0;
    float mass = 0.0;
    for (unsigned i = 0; i < nAtoms; ++i)
    {
        ar >> pos >> mass;
        domain.addAtom(Atom(pos, mass)); // this will increment mNumAtoms
    }
    return ar;
}