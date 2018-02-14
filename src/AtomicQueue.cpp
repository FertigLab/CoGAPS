#include "AtomicQueue.h"

// TODO have private random number generator for consistency with benchmark

// consider having separate class that generates the proposal sequence
// and selects the atoms, but does not calculate the mass .... separate
// class with matrix access computes the mass - also allows for entire
// proposal specific code to be together (no more conditionals on nChanges)

// acceptance of proposals needs it's own unit as well - considering
// how much time is spent updating matrices - the updates would be independent
// of any currently running calculations by design. Find a way to lock
// the individual vectors of the matrix without locking the entire matrix.

// consider changing the data layout to have the row matrices (D,AP,S) near
// the P matrix and the col matrices (D,AP,S) near the A matrix

// consider allocating a single continuous chunk for each matric instead
// of individual vectors - does this mess up locking?

// look into delaying the A/P matrix update until after the full update step
// - only the AP updates are neccesary

// no need to separate alphaParameters and deltaLL into separate workers - 
// they are serial and touch much of the same memory

// need to ensure proposals dont effect each other
// BIRTH - mark off locations of birth
// DEATH - can't free up death locations, neglible effect
// MOVE - must move within neighbor boundary, if a neighbor could be moved
//  or deleted then this is held up - also know where birth locations are
// EXCHANGE - selects right most atom - uneffected by move, by death of this
//  atom or birth of an atom in between can effect

// Neighbors
// DEATH - merge neighbors
// MOVE/RESIZE - nothing changes
// EXCHANGE - nothing changes
// BIRTH - need to find left or right neighbor
// store neighbor information, then everything constant besides insertion
// store in map makes insertion logN

AtomicQueue::AtomicQueue(char label, unsigned nrow, unsigned ncol)
    :
mLabel(label), mNumAtoms(0),
mMaxNumAtoms(std::numeric_limits<uint64_t>::max()),
mNumRows(nrow), mNumCols(ncol), mNumBins(nrow * ncol),
mBinSize(std::numeric_limits<uint64_t>::max() / (nrow * ncol)),
mAlpha(alpha), mLambda(lambda)
{}

// O(logN)
void AtomicQueue::addAtom(Atom atom)
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
        mMaxAtoms--; // keep queue up to date
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

static isInVector(std::vector<unsigned> &vec, unsigned n)
{
    return std::find(vec.begin(), vec.end(), n) != vec.end();
}

bool AtomicSupport::addProposal() const
{
    if (minAtoms == 0 && maxAtoms > 0 || minAtoms < 2 && maxAtoms >= 2)
    {
        return false;
    }

    AtomicProposal proposal;
    if (maxAtoms == 0)
    {
        proposal = proposeBirth();
        minAtoms++;
        maxAtoms++;
    }
    else
    {
        float bdProb = maxAtoms < 2 ? 0.6667f : 0.5f;
        float u1 = gaps::random::uniform(); // cache these values if a failure
        float u2 = gaps::random::uniform(); // happens, needed to prevent bias
        float lowerBound = deleteProb(minAtoms);
        float upperBound = deleteProb(maxAtoms);
        if (u1 < bdProb && u2 > upperBound)
        {
            proposal = proposeBirth();
            minAtoms++;
            maxAtoms++;
        }
        else if (u1 < bdProb && u2 < lowerBound)
        {
            proposal = proposeDeath();
            minAtoms--;
        }
        else if (u1 >= bdProb && (maxAtoms < 2 || u1 < 0.75f))
        {
            proposal = proposeMove();
        }
        else if (u1 >= bdProb)
        {
            proposal = proposeExchange();
            minAtoms--;
        }
        else
        {
            return false;
        }
    }
    
    // used rows should be an hash table so that find is O(1)
    // can allocate is one chunk of memory since nRows/nCols is known
    if (isInVector(usedRows, getRow(proposal.pos1)))
    {
        return false;
    }
    if (proposal.nChanges > 1)
    {
        if (isInVector(usedRows, getRow(proposal.pos2)))
        {
            return false;
        }
        usedRows.push_back(getRow(proposal.pos2));
    }
    usedRows.push_back(getRow(proposal.pos1));
    mProposalQueue.push(proposal);
    return true;
}

void AtomicSupport::populateQueue(unsigned limit)
{
    unsigned nIter = 0;
    while (nIter++ < limit && addProposal());
}

AtomicProposal AtomicSupport::popQueue()
{
    AtomicProposal prop = mProposalQueue.front();
    mProposalQueue.pop();
    return prop;
}

MatrixChange AtomicSupport::acceptProposal(const AtomicProposal &prop,
MatrixChange ch)
{
    ch.delta1 = updateAtomMass(prop.pos1, prop.delta1);
    ch.delta2 = (prop.nChanges > 1) ? updateAtomMass(prop.pos2, prop.delta2)
        : ch.delta2;
    return ch;
}

void AtomicSupport::rejectProposal(const AtomicProposal &prop)
{
    if (prop.label == 'D' || prop.label == 'E')
    {
        mMinAtoms++;
    }
}

unsigned AtomicSupport::size() const
{
    return mProposalQueue.size();
}

bool AtomicSupport::empty() const
{
    return mProposalQueue.empty();
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