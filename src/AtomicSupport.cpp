#include "AtomicSupport.h"

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <vector>
#include <limits>

#include <sstream>
#include <string>
#include <istream>
#include <iterator>
#include <fstream>

#define EPSILON 1E-10

AtomicSupport::AtomicSupport(uint64_t nrow, uint64_t ncol)
    : mNumRows(nrow), mNumCols(ncol), mNumBins(nrow * ncol),
      mNumAtoms(0), mTotalMass(0.0),
      mMaxNumAtoms(std::numeric_limits<uint64_t>::max()),
      mBinSize(std::numeric_limits<uint64_t>::max() / (nrow * ncol))
{}

uint64_t AtomicSupport::getRow(uint64_t pos) const
{
    uint64_t binNum = std::max(pos / mBinSize, mNumBins - 1);
    return binNum / mNumCols;
}

uint64_t AtomicSupport::getCol(uint64_t pos) const
{
    uint64_t binNum = std::max(pos / mBinSize, mNumBins - 1);
    return binNum % mNumCols;
}

uint64_t AtomicSupport::randomFreePosition() const
{
    uint64_t pos = 0;
    do
    {
        pos = Random::uniform64();
    } while (mAtomicDomain.count(pos));
    return pos;
}

uint64_t AtomicSupport::randomAtomPosition() const
{
    uint64_t num = Random::uniform64(0, mNumAtoms - 1);
    std::map<uint64_t, double>::const_iterator it = mAtomicDomain.begin();
    std::advance(it, num);
    return it->first;
}

AtomicProposal AtomicSupport::proposeDeath() const
{
    uint64_t location = randomAtomPosition();
    uint64_t mass = mAtomicDomain.at(location);
    return AtomicProposal('D', location, -mass);
}

AtomicProposal AtomicSupport::proposeBirth() const
{
    uint64_t location = randomFreePosition();
    uint64_t mass = Random::exponential(mLambda);
    return AtomicProposal('B', location, mass);
}

// move atom between adjacent atoms
AtomicProposal AtomicSupport::proposeMove() const
{
    // get random atom
    uint64_t location = randomAtomPosition();
    uint64_t mass = mAtomicDomain.at(location);

    // get an iterator to this atom
    std::map<uint64_t, double>::const_iterator it;
    it = mAtomicDomain.find(location);
    uint64_t left = (it == mAtomicDomain.begin() ? it->first : (++it)->first);
    if (++it == mAtomicDomain.end())
        --it;
    uint64_t right = it->first;

    uint64_t newLocation = Random::uniform64(left, right);

    return AtomicProposal('M', location, -mass, newLocation, mass);
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
    double pupper = Random::p_gamma(mass1 + mass2, 2.0, 1.0 / mLambda);
    double newMass = Random::q_gamma(Random::uniform(0.0, pupper), 2.0,
        1.0 / mLambda);

    // preserve total mass
    return mass1 > mass2 ?
        AtomicProposal('E', pos1, newMass - mass1, pos2, mass1 - newMass)
        : AtomicProposal('E', pos1, mass2 - newMass, pos2, newMass - mass2);
}

void AtomicSupport::updateAtomMass(uint64_t pos, double delta)
{
    if (mAtomicDomain.count(pos)) // update atom if it exists
    {
        mAtomicDomain.at(pos) += delta;
        mTotalMass += delta;
    }
    else // create a new atom if it doesn't exist
    {
        mAtomicDomain.insert(std::pair<uint64_t, double>(pos, delta));
        mTotalMass += delta;
        mNumAtoms++;
    }
    if (mAtomicDomain.at(pos) <= EPSILON) // check if atom is too small
    {
        mTotalMass -= mAtomicDomain.at(pos);
        mAtomicDomain.erase(pos);
        mNumAtoms--;
    }
}

#define GEOM_F(x) (((double) (x)) / ((double) (x + 1)))
#define POISS_F(x,y) ((double)((x) * (y)) / (double)((x) - (y)))
AtomicProposal AtomicSupport::makeProposal() const
{
    if (mNumAtoms == 0)
    {
        return proposeBirth();
    }

    double unif = Random::uniform();
    if (mNumAtoms < 2 && unif <= 0.6667 || unif <= 0.5) // birth/death
    {
        if (mNumAtoms >= mMaxNumAtoms)
        {
            return proposeDeath();
        }

        double pDelete = 0.5; // default uniform prior for mAlpha == 0
        if (mAlpha > 0) // poisson prior
        {
            pDelete = 1 + POISS_F(mMaxNumAtoms, mNumAtoms) / (mAlpha * mNumBins);
        }
        else if (mAlpha < 0) // geometric prior
        {
            pDelete = 1 + GEOM_F(mNumAtoms) / GEOM_F(-mAlpha * mNumBins);
        }
        return Random::uniform() < pDelete ? proposeDeath() : proposeBirth();
    }
    return mNumAtoms < 2 || unif <= 0.75 ? proposeMove() : proposeExchange();
}

void AtomicSupport::acceptProposal(const AtomicProposal &prop)
{
    updateAtomMass(prop.pos1, prop.deltaMass1);
    if (prop.nChanges > 1)
    {
        updateAtomMass(prop.pos2, prop.deltaMass2);
    }
}

std::vector<ElementChange> AtomicSupport::getElementChange(const AtomicProposal &prop) const
{
    std::vector<ElementChange> vec;
    vec.push_back(ElementChange(getRow(prop.pos1), getCol(prop.pos1), prop.deltaMass1));

    if (prop.nChanges > 1)
    {
        vec.push_back(ElementChange(getRow(prop.pos2), getCol(prop.pos2), prop.deltaMass2));
    }
    return vec;
}

void AtomicSupport::write(const std::string &path, bool append) const
{
    std::ios_base::openmode mode = append ? std::ios::out | std::ios::app
        : std::ios::out;
    std::ofstream outputFile(path.c_str(), mode);

    std::map<uint64_t, double>::const_iterator iter;
    for (iter = mAtomicDomain.begin(); iter != mAtomicDomain.end(); ++iter)
    {
        outputFile << setiosflags(std::ios::right) << std::setw(25) << iter->first
            << ' ' << std::setw(15) << iter->second << '\n';
    }
    outputFile.close();
}