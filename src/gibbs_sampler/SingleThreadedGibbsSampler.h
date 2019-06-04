#ifndef __COGAPS_SINGLE_THREADED_GIBBS_SAMPLER_H__
#define __COGAPS_SINGLE_THREADED_GIBBS_SAMPLER_H__

#include "AlphaParameters.h"
#include "../atomic/AtomicDomain.h"
#include "../data_structures/Matrix.h"
#include "../math/Math.h"
#include "../math/VectorMath.h"
#include "../math/MatrixMath.h"
#include "../math/Random.h"
#include "../GapsParameters.h"

#include <vector>
#include <limits>
#include <cmath>
#include <stdint.h>

////////////////////////// SingleThreadedGibbsSampler Interface //////////////////////////

class GapsStatistics;

template <class DataModel>
class SingleThreadedGibbsSampler;

template <class DataModel>
Archive& operator<<(Archive &ar, const SingleThreadedGibbsSampler<DataModel> &sampler);

template <class DataModel>
Archive& operator>>(Archive &ar, SingleThreadedGibbsSampler<DataModel> &sampler);

template <class DataModel>
class SingleThreadedGibbsSampler : public DataModel
{
public:

    template <class DataType>
    SingleThreadedGibbsSampler(const DataType &data, bool transpose, bool subsetRows,
        float alpha, float maxGibbsMass, const GapsParameters &params,
        GapsRandomState *randState);

    unsigned nAtoms() const;
    float getAverageQueueLength() const;

    void update(unsigned nSteps, unsigned nThreads);

    friend Archive& operator<< <DataModel> (Archive &ar, const SingleThreadedGibbsSampler &s);
    friend Archive& operator>> <DataModel> (Archive &ar, SingleThreadedGibbsSampler &s);

private:

    AtomicDomain mDomain; // data structure providing access to atoms

    mutable GapsRng mRng;

    uint64_t mNumBins;
    uint64_t mBinLength;
    uint64_t mNumPatterns;

    double mDomainLength; // length of entire atomic domain
    double mAlpha;

    char getUpdateType() const;
    void birth();
    void death();
    void move();
    void exchange();
};

/////////////////////// SingleThreadedGibbsSampler Implementation ////////////////////////

template <class DataModel>
template <class DataType>
SingleThreadedGibbsSampler<DataModel>::SingleThreadedGibbsSampler(const DataType &data,
bool transpose, bool subsetRows, float alpha, float maxGibbsMass,
const GapsParameters &params, GapsRandomState *randState)
    :
DataModel(data, transpose, subsetRows, params, alpha),
mDomain(DataModel::nElements(), randState),
mRng(randState),
mNumBins(DataModel::nElements()),
mBinLength(std::numeric_limits<uint64_t>::max() / (DataModel::nElements())),
mNumPatterns(DataModel::nPatterns()),
mDomainLength(mBinLength * DataModel::nElements()),
mAlpha(alpha)
{}

template <class DataModel>
unsigned SingleThreadedGibbsSampler<DataModel>::nAtoms() const
{
    return mDomain.size();
}

template <class DataModel>
float SingleThreadedGibbsSampler<DataModel>::getAverageQueueLength() const
{
    return 0.f;
}

template <class DataModel>
char SingleThreadedGibbsSampler<DataModel>::getUpdateType() const
{
    if (mDomain.size() < 2)
    {
        return 'B'; // always birth when no atoms exist
    }

    float u1 = mRng.uniform();
    if (u1 < 0.5f)
    {
        double nAtoms = static_cast<double>(mDomain.size());
        double numer = nAtoms * mDomainLength;
        float deathProb = numer / (numer + mAlpha * mNumBins * (mDomainLength - nAtoms));
        return mRng.uniform() < deathProb ? 'D' : 'B';
    }
    return u1 < 0.75f ? 'M' : 'E';
}

template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::update(unsigned nSteps, unsigned nThreads)
{
    for (unsigned i = 0; i < nSteps; ++i)
    {
        switch (getUpdateType())
        {
            case 'B': birth();    break;
            case 'D': death();    break;
            case 'M': move();     break;
            case 'E': exchange(); break;
        }
    }

    GAPS_ASSERT(mDomain.isSorted());
}

// add an atom at a random position, calculate mass either with an
// exponential distribution or with the gibbs mass distribution
template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::birth()
{
    // get random open position in atomic domain, calculate row and col of the position
    uint64_t pos = mDomain.randomFreePosition();
    unsigned row = (pos / mBinLength) / mNumPatterns;
    unsigned col = (pos / mBinLength) % mNumPatterns;

    // try to get mass using gibbs, resort to exponential if needed
    OptionalFloat mass = DataModel::canUseGibbs(col)
        ? DataModel::sampleBirth(row, col, &mRng)
        : mRng.exponential(DataModel::lambda());

    // accept mass as long as gibbs succeded or it's non-zero
    if (mass.hasValue() && mass.value() > gaps::epsilon)
    {
        mDomain.insert(pos, mass.value());
        DataModel::changeMatrix(row, col, mass.value());
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass distribution
template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::death()
{
    // select atom at random and calculate it's row and col
    Atom *atom = mDomain.randomAtom();
    unsigned row = (atom->pos / mBinLength) / mNumPatterns;
    unsigned col = (atom->pos / mBinLength) % mNumPatterns;

    // try to do a rebirth in the place of this atom
    if (DataModel::canUseGibbs(col))
    {
        OptionalFloat mass = DataModel::sampleDeathAndRebirth(row, col,
            -1.f * atom->mass, &mRng);
        if (mass.hasValue() && mass.value() >= gaps::epsilon)
        {
            DataModel::safelyChangeMatrix(row, col, mass.value() - atom->mass);
            atom->mass = mass.value();
            return;
        }
    }

    // if rebirth fails, then kill off atom
    DataModel::safelyChangeMatrix(row, col, -1.f * atom->mass);
    mDomain.erase(atom->pos);
}

// move mass from src to dest in the atomic domain
template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::move()
{
    // select atom at random and get it's right and left neighbors
    AtomNeighborhood hood = mDomain.randomAtomWithNeighbors();
    uint64_t lbound = hood.hasLeft() ? hood.left->pos : 0;
    uint64_t rbound = hood.hasRight() ? hood.right->pos :
        static_cast<uint64_t>(mDomainLength);

    // randomly select new position to move to and calculation the row and col of it
    uint64_t pos = mRng.uniform64(lbound + 1, rbound - 1);
    Atom *atom = hood.center;
    unsigned r1 = (atom->pos / mBinLength) / mNumPatterns;
    unsigned c1 = (atom->pos / mBinLength) % mNumPatterns;
    unsigned r2 = (pos / mBinLength) / mNumPatterns;
    unsigned c2 = (pos / mBinLength) % mNumPatterns;

    // automatically accept move if it keeps atom in the same matrix element
    if (r1 == r2 && c1 == c2)
    {
        atom->pos = pos;
        return;
    }
    
    // conditionally accept move based on change to likelihood
    float deltaLL = DataModel::deltaLogLikelihood(r1, c1, r2, c2, atom->mass);
    if (std::log(mRng.uniform()) < deltaLL)
    {
        atom->pos = pos;
        DataModel::safelyChangeMatrix(r1, c1, -atom->mass);
        DataModel::changeMatrix(r2, c2, atom->mass);
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small
template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::exchange()
{
    // select atom at random and get it's right neighbor
    AtomNeighborhood hood = mDomain.randomAtomWithRightNeighbor();
    Atom *atom1 = hood.center;
    Atom *atom2 = hood.hasRight() ? hood.right : mDomain.front();

    // calculate row and col of the atom and it's neighbor
    unsigned r1 = (atom1->pos / mBinLength) / mNumPatterns;
    unsigned c1 = (atom1->pos / mBinLength) % mNumPatterns;
    unsigned r2 = (atom2->pos / mBinLength) / mNumPatterns;
    unsigned c2 = (atom2->pos / mBinLength) % mNumPatterns;

    // ignore exhanges in the same bin
    if ((r1 != r2 || c1 != c2) && DataModel::canUseGibbs(c1, c2))
    {
        OptionalFloat mass = DataModel::sampleExchange(r1, c1, atom1->mass,
            r2, c2, atom2->mass, &mRng);
        float newMass1 = atom1->mass + mass.value();
        float newMass2 = atom2->mass - mass.value();
        if (mass.hasValue() && newMass1 > gaps::epsilon && newMass2 > gaps::epsilon)
        {
            DataModel::safelyChangeMatrix(r1, c1, newMass1 - atom1->mass);
            DataModel::safelyChangeMatrix(r2, c2, newMass2 - atom2->mass);
            atom1->mass = newMass1;
            atom2->mass = newMass2;
            return;
        }
    }
}

template <class DataModel>
Archive& operator<<(Archive &ar, const SingleThreadedGibbsSampler<DataModel> &s)
{
    operator<<(ar, static_cast<const DataModel&>(s)) << s.mDomain << s.mNumBins
        << s.mBinLength << s.mNumPatterns << s.mDomainLength << s.mAlpha;
    return ar;
}

template <class DataModel>
Archive& operator>>(Archive &ar, SingleThreadedGibbsSampler<DataModel> &s)
{
    operator>>(ar, static_cast<DataModel&>(s)) << s.mDomain << s.mNumBins
        << s.mBinLength << s.mNumPatterns << s.mDomainLength << s.mAlpha;
    return ar;
}

#endif // __COGAPS_SINGLE_THREADED_GIBBS_SAMPLER_H__
