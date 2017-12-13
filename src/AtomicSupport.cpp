#include "AtomicSupport.h"

static const double EPSILON = 1.e-10;

AtomicSupport::AtomicSupport(char label, uint64_t nrow, uint64_t ncol,
double alpha, double lambda)
    :
mNumRows(nrow), mNumCols(ncol), mNumBins(nrow * ncol),
mNumAtoms(0), mTotalMass(0.0), mLabel(label), mAlpha(alpha), 
mMaxNumAtoms(std::numeric_limits<uint64_t>::max()), mLambda(lambda),
mBinSize(std::numeric_limits<uint64_t>::max() / (nrow * ncol))
{}

uint64_t AtomicSupport::getRow(uint64_t pos) const
{
    uint64_t binNum = std::min(pos / mBinSize, mNumBins - 1);
    return binNum / mNumCols;
}

uint64_t AtomicSupport::getCol(uint64_t pos) const
{
    uint64_t binNum = std::min(pos / mBinSize, mNumBins - 1);
    return binNum % mNumCols;
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
    std::map<uint64_t, double>::const_iterator it = mAtomicDomain.begin();
    std::advance(it, num);
    return it->first;
}

AtomicProposal AtomicSupport::proposeDeath() const
{
    uint64_t location = randomAtomPosition();
    double mass = mAtomicDomain.at(location);
    return AtomicProposal(mLabel, 'D', location, -mass);
}

AtomicProposal AtomicSupport::proposeBirth() const
{
    uint64_t location = randomFreePosition();
    double mass = gaps::random::exponential(mLambda);
    return AtomicProposal(mLabel, 'B', location, std::max(mass, EPSILON));
}

// move atom between adjacent atoms
AtomicProposal AtomicSupport::proposeMove() const
{
    // get random atom
    uint64_t location = randomAtomPosition();
    double mass = mAtomicDomain.at(location);

    // get an iterator to this atom
    std::map<uint64_t, double>::const_iterator it;
    it = mAtomicDomain.find(location);
    uint64_t left = (it == mAtomicDomain.begin() ? it->first : (++it)->first);
    if (++it == mAtomicDomain.end())
        --it;
    uint64_t right = it->first;

    uint64_t newLocation = gaps::random::uniform64(left, right);

    return AtomicProposal(mLabel, 'M', location, -mass, newLocation, mass);
}

// exchange with adjacent atom to the right
AtomicProposal AtomicSupport::proposeExchange() const
{
    // get random atom
    uint64_t pos1 = randomAtomPosition();
    double mass1 = mAtomicDomain.at(pos1);

    // find atom to the right, wrap around the end
    std::map<uint64_t, double>::const_iterator it;
    it = ++(mAtomicDomain.find(pos1));
    if (it == mAtomicDomain.end())
    {
        it = mAtomicDomain.begin();
    }

    // get adjacent atom info
    uint64_t pos2 = it->first;
    double mass2 = it->second;

    // calculate new mass
    double pupper = gaps::random::p_gamma(mass1 + mass2, 2.0, 1.0 / mLambda);
    double newMass = gaps::random::q_gamma(gaps::random::uniform(0.0, pupper),
        2.0, 1.0 / mLambda);

    // preserve total mass
    return mass1 > mass2 ?
        AtomicProposal(mLabel, 'E', pos1, newMass - mass1, pos2, mass1 - newMass)
        : AtomicProposal(mLabel, 'E', pos1, mass2 - newMass, pos2, newMass - mass2);
}

double AtomicSupport::updateAtomMass(uint64_t pos, double delta)
{
    if (mAtomicDomain.count(pos)) // update atom if it exists
    {
        mAtomicDomain.at(pos) += delta;
    }
    else // create a new atom if it doesn't exist
    {
        mAtomicDomain.insert(std::pair<uint64_t, double>(pos, delta));
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

    double unif = gaps::random::uniform();
    if ((mNumAtoms < 2 && unif <= 0.6667) || unif <= 0.5) // birth/death
    {
        if (mNumAtoms >= mMaxNumAtoms)
        {
            return proposeDeath();
        }

        double pDelete = 0.5; // default uniform prior for mAlpha == 0
        if (mAlpha > 0) // poisson prior
        {
            double dMax = (double)mMaxNumAtoms;
            double dNum = (double)mNumAtoms;
            double maxTerm = (dMax - dNum) / dMax;

            pDelete = dNum / (dNum + mAlpha * (double)mNumBins * maxTerm);
        }
        else if (mAlpha < 0) // geometric prior
        {
            double c = -mAlpha * mNumBins / (-mAlpha * mNumBins + 1.0);
            pDelete = (double)mNumAtoms / ((mNumAtoms + 1) * c + (double)mNumAtoms);
        }
        return gaps::random::uniform() < pDelete ? proposeDeath() : proposeBirth();
    }
    return (mNumAtoms < 2 || unif >= 0.75) ? proposeMove() : proposeExchange();
}

MatrixChange AtomicSupport::acceptProposal(const AtomicProposal &prop)
{
#ifdef GAPS_DEBUG
    if (prop.label != mLabel)
    {
        throw std::runtime_error("domain mismatch");
    }
#endif

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
#ifdef GAPS_DEBUG
    if (prop.label != mLabel)
    {
        throw std::runtime_error("domain mismatch");
    }
#endif

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