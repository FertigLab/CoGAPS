#include "AtomicSupport.h"

static const float EPSILON = 1.e-10;

AtomicSupport::AtomicSupport(char label, uint64_t nrow, uint64_t ncol,
float alpha, float lambda)
    :
mLabel(label), mNumAtoms(0),
mMaxNumAtoms(std::numeric_limits<uint64_t>::max()),
mTotalMass(0.0), mNumRows(nrow), mNumCols(ncol), mNumBins(nrow * ncol),
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

uint64_t AtomicSupport::getRow(uint64_t pos) const
{
    uint64_t binNum = std::min(pos / mBinSize, mNumBins - 1);
    return mLabel == 'A' ? binNum / mNumCols : binNum % mNumRows;
}

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

AtomicProposal AtomicSupport::proposeDeath() const
{
    uint64_t location = randomAtomPosition();
    float mass = mAtoms[mAtomPositions.at(location)].mass;
    return AtomicProposal(mLabel, 'D', location, -mass);
}

AtomicProposal AtomicSupport::proposeBirth() const
{
    uint64_t location = randomFreePosition();
    float mass = gaps::random::exponential(mLambda);
    return AtomicProposal(mLabel, 'B', location, mass);
}

// move atom between adjacent atoms
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
AtomicProposal AtomicSupport::proposeExchange() const
{
    uint64_t pos1 = randomAtomPosition();
    AtomNeighbors nbs = getNeighbors(pos1);
    float mass1 = at(pos1);

    // find right neighbor, wrap around to beginning if this atom is end
    uint64_t pos2 = nbs.right.pos != pos1 ? nbs.right.pos : mAtomPositions.begin()->first;
    float mass2 = at(pos2);

    // calculate new mass
    float pupper = gaps::random::p_gamma(mass1 + mass2, 2.0, 1.0 / mLambda);
    float newMass = gaps::random::q_gamma(gaps::random::uniform(0.0, pupper),
        2.0, 1.0 / mLambda);

    // calculate mass changes
    float delta1 = mass1 > mass2 ? newMass - mass1 : mass2 - newMass;
    float delta2 = mass2 > mass1 ? newMass - mass2 : mass1 - newMass;

    // preserve total mass
    return AtomicProposal(mLabel, 'E', pos1, delta1, pos2, delta2);
}

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

    if (mAtoms[mAtomPositions.at(pos)].mass < EPSILON) // check if atom is too small
    {
        delta -= at(pos);
        removeAtom(pos);
    }
    mTotalMass += delta;
    return delta;
}

AtomicProposal AtomicSupport::makeProposal() const
{
    if (mNumAtoms == 0)
    {
        return proposeBirth();
    }

    float unif = gaps::random::uniform();
    if ((mNumAtoms < 2 && unif <= 0.6667f) || unif <= 0.5f) // birth/death
    {
        if (mNumAtoms >= mMaxNumAtoms)
        {
            return proposeDeath();
        }

        float pDelete = 0.5; // default uniform prior for mAlpha == 0
        if (mAlpha > 0) // poisson prior
        {
            float dMax = (float)mMaxNumAtoms;
            float dNum = (float)mNumAtoms;
            float maxTerm = (dMax - dNum) / dMax;

            pDelete = dNum / (dNum + mAlpha * (float)mNumBins * maxTerm);
        }
        else if (mAlpha < 0) // geometric prior
        {
            float c = -mAlpha * mNumBins / (-mAlpha * mNumBins + 1.0);
            pDelete = (float)mNumAtoms / ((mNumAtoms + 1) * c + (float)mNumAtoms);
        }
        return gaps::random::uniform() < pDelete ? proposeDeath() : proposeBirth();
    }
    return (mNumAtoms < 2 || unif < 0.75) ? proposeMove() : proposeExchange();
}

MatrixChange AtomicSupport::acceptProposal(const AtomicProposal &prop)
{
    MatrixChange change = getMatrixChange(prop);
    change.delta1 = updateAtomMass(prop.pos1, prop.delta1);
    if (prop.nChanges > 1)
    {
        change.delta2 = updateAtomMass(prop.pos2, prop.delta2);
    }
    return change;
}

MatrixChange AtomicSupport::getMatrixChange(const AtomicProposal &prop) const
{
    if (prop.nChanges > 1)
    {
        return MatrixChange(prop.label, getRow(prop.pos1), getCol(prop.pos1), prop.delta1,
            getRow(prop.pos2), getCol(prop.pos2), prop.delta2);
    }
    else
    {   
        return MatrixChange(prop.label, getRow(prop.pos1), getCol(prop.pos1), prop.delta1);
    }
}

Archive& operator<<(Archive &ar, AtomicSupport &domain)
{
    ar << domain.mLabel << domain.mNumAtoms << domain.mMaxNumAtoms
        << domain.mTotalMass << domain.mNumRows << domain.mNumCols
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
    ar >> domain.mLabel >> nAtoms >> domain.mMaxNumAtoms >> domain.mTotalMass
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