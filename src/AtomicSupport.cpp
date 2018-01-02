#include "AtomicSupport.h"

static const float EPSILON = 1.e-10;

AtomicSupport::AtomicSupport(char label, uint64_t nrow, uint64_t ncol,
float alpha, float lambda)
    :
mNumRows(nrow), mNumCols(ncol), mNumBins(nrow * ncol),
mNumAtoms(0), mTotalMass(0.0), mLabel(label), mAlpha(alpha), 
mMaxNumAtoms(std::numeric_limits<uint64_t>::max()), mLambda(lambda),
mBinSize(std::numeric_limits<uint64_t>::max() / (nrow * ncol))
{}

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

uint64_t AtomicSupport::randomFreePosition() const
{
    uint64_t pos = 0;
    do
    {
        pos = gaps::random::uniform64();
    } while (mAtomicDomain.count(pos) > 0);
    return pos;
}

uint64_t AtomicSupport::randomAtomPosition() const
{
    uint64_t num = gaps::random::uniform64(0, mNumAtoms - 1);
    std::map<uint64_t, float>::const_iterator it = mAtomicDomain.begin();
    std::advance(it, num);
    return it->first;
}

AtomicProposal AtomicSupport::proposeDeath() const
{
    uint64_t location = randomAtomPosition();
    float mass = mAtomicDomain.at(location);
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
    // get random atom
    uint64_t location = randomAtomPosition();
    float mass = mAtomicDomain.at(location);

    // get an iterator to this atom
    std::map<uint64_t, float>::const_iterator it, left, right;
    it = mAtomicDomain.find(location);
    left = it; right = it;
    uint64_t rbound = mMaxNumAtoms - 1, lbound = 0;
    if (left != mAtomicDomain.begin())
    {
        --left;
        lbound = left->first;
    }
    if (++right != mAtomicDomain.end())
    {
        rbound = right->first;
    }

    uint64_t newLocation = gaps::random::uniform64(lbound, rbound);
    return AtomicProposal(mLabel, 'M', location, -mass, newLocation, mass);
}

// exchange with adjacent atom to the right
AtomicProposal AtomicSupport::proposeExchange() const
{
    // get random atom
    uint64_t pos1 = randomAtomPosition();
    float mass1 = mAtomicDomain.at(pos1);

    // find atom to the right, wrap around the end
    std::map<uint64_t, float>::const_iterator it;
    it = ++(mAtomicDomain.find(pos1));
    if (it == mAtomicDomain.end())
    {
        it = mAtomicDomain.begin();
    }

    // get adjacent atom info
    uint64_t pos2 = it->first;
    float mass2 = it->second;

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

float AtomicSupport::updateAtomMass(char type, uint64_t pos, float delta)
{
    if (mAtomicDomain.count(pos)) // update atom if it exists
    {
        mAtomicDomain.at(pos) += delta;
    }
    else // create a new atom if it doesn't exist
    {
        mAtomicDomain.insert(std::pair<uint64_t, float>(pos, delta));
        mNumAtoms++;
    }

    if (mAtomicDomain.at(pos) < EPSILON) // check if atom is too small
    {
        delta -= mAtomicDomain.at(pos);
        mAtomicDomain.erase(pos);
        mNumAtoms--;
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
    if ((mNumAtoms < 2 && unif <= 0.6667) || unif <= 0.5) // birth/death
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
    change.delta1 = updateAtomMass(prop.type, prop.pos1, prop.delta1);
    if (prop.nChanges > 1)
    {
        change.delta2 = updateAtomMass(prop.type, prop.pos2, prop.delta2);
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