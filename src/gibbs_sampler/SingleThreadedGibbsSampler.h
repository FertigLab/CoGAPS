#ifndef __COGAPS_SINGLE_THREADED_GIBBS_SAMPLER_H__
#define __COGAPS_SINGLE_THREADED_GIBBS_SAMPLER_H__

#include "AlphaParameters.h"
#include "../data_structures/Matrix.h"
#include "../math/Math.h"
#include "../math/VectorMath.h"
#include "../math/MatrixMath.h"
#include "../math/Random.h"
#include "../GapsParameters.h"

#include <cmath>
#include <limits>
#include <stdint.h>
#include <vector>

#include "../atomic/AtomicDomain.h"
typedef Atom AtomType;
typedef AtomNeighborhood AtomNeighborhoodType;
typedef AtomicDomain AtomicDomainType;

////////////////////////// SingleThreadedGibbsSampler Interface //////////////////////////

class GapsStatistics;

template <class DataModel>
class SingleThreadedGibbsSampler;

template <class DataModel>
Archive& operator<<(Archive &ar, const SingleThreadedGibbsSampler<DataModel> &s);

template <class DataModel>
Archive& operator>>(Archive &ar, SingleThreadedGibbsSampler<DataModel> &s);

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
    char getUpdateType() const;
    void birth();
    void death();
    void move();
    void exchange();

    AtomicDomainType mDomain; // data structure providing access to atoms
    mutable GapsRng mRng;
    uint64_t mNumBins;
    uint64_t mBinLength;
    uint64_t mNumPatterns;
    double mDomainLength; // length of entire atomic domain
    double mAlpha;
};

/////////////////////// SingleThreadedGibbsSampler Implementation ////////////////////////

template <class DataModel>
template <class DataType>
SingleThreadedGibbsSampler<DataModel>::SingleThreadedGibbsSampler(const DataType &data,
bool transpose, bool subsetRows, float alpha, float maxGibbsMass,
const GapsParameters &params, GapsRandomState *randState)
    :
DataModel(data, transpose, subsetRows, params, alpha, maxGibbsMass),
mDomain(DataModel::nElements()),
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
    // Move and exchange require existing atoms and neighboring atom positions.
    // With fewer than two atoms, force birth so early model-fitting states do
    // not create degenerate proposals with undefined conditional distributions.
    if (mDomain.size() < 2)
    {
        return 'B'; // always birth when no atoms exist
    }

    // Choose among the four reversible-jump moves used by the atomic prior.
    // Birth/death change the number of atoms; move/exchange preserve it.
    float u1 = mRng.uniform();
    if (u1 < 0.5f)
    {
        double nAtoms = static_cast<double>(mDomain.size());
        double numer = nAtoms * mDomainLength;
        float deathProb = numer / (numer + mAlpha * mNumBins * (mDomainLength - nAtoms));
        // Within the birth/death branch, choose death with probability implied
        // by the atom-count prior; otherwise propose birth.
        return mRng.uniform() < deathProb ? 'D' : 'B';
    }
    // Split the remaining proposal interval between move and exchange.
    return u1 < 0.75f ? 'M' : 'E';
}

template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::update(unsigned nSteps, unsigned nThreads) // NOLINT
{
    // Sequential sampler path: each proposal is generated, evaluated, and
    // applied before the next proposal is drawn.
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
}

// Birth adds an atom at a random free position in the atomic domain. Early in
// model fitting, the conditional Gibbs distribution may be undefined because
// its variance term would require division by zero. In those cases, the atom
// mass is sampled from the exponential prior so the full proposal can still be
// evaluated.
template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::birth()
{
    // Map the atomic position to the matrix row and pattern column it updates.
    uint64_t pos = mDomain.randomFreePosition(&mRng);
    unsigned row = (pos / mBinLength) / mNumPatterns;
    unsigned col = (pos / mBinLength) % mNumPatterns;

    // Sample mass using Gibbs when the conditional distribution is defined;
    // fall back to the prior when the conditional variance would be degenerate.
    OptionalFloat mass = DataModel::canUseGibbs(col)
        ? DataModel::sampleBirth(row, col, &mRng)
        : mRng.exponential(DataModel::lambda());

    // Accept the birth when the proposed mass is valid and nonzero.
    if (mass.hasValue() && mass.value() > gaps::epsilon)
    {
        mDomain.insert(pos, mass.value());
        DataModel::changeMatrix(row, col, mass.value());
    }
}

// Death removes an atom by first subtracting its mass from the matrix, then
// attempting a rebirth at the same matrix element. When the conditional Gibbs
// distribution is defined, it proposes the rebirth mass; otherwise the original
// mass is retained so the full death/rebirth proposal can still be evaluated.
template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::death()
{
    // Select an atom and map its position to the matrix row and pattern column.
    AtomType *atom = mDomain.randomAtom(&mRng);
    unsigned row = (atom->pos() / mBinLength) / mNumPatterns;
    unsigned col = (atom->pos() / mBinLength) % mNumPatterns;
    // Determine the mass to attempt rebirth with.
    float rebirthMass = atom->mass(); // default rebirth mass == no change to atom
    AlphaParameters alpha = DataModel::alphaParametersWithChange(row, col, -1.f * atom->mass())
        * DataModel::annealingTemp();
    if (DataModel::canUseGibbs(col))
    {
        OptionalFloat gMass = gibbsMass(alpha, 0.f, DataModel::maxGibbsMass(), &mRng,
            DataModel::lambda());
        if (gMass.hasValue())
        {
            rebirthMass = gMass.value();
        }
    }
    // Handle accept/reject of the rebirth step.
    float deltaLL = rebirthMass * (alpha.s_mu - alpha.s * rebirthMass / 2.f);
    if (std::log(mRng.uniform()) < deltaLL) // accept
    {
        if (rebirthMass != atom->mass())
        {
            DataModel::safelyChangeMatrix(row, col, rebirthMass - atom->mass());
            atom->updateMass(rebirthMass);
        }
    }
    else // reject
    {
        DataModel::safelyChangeMatrix(row, col, -1.f * atom->mass());
        mDomain.erase(atom);
    }
}

// Move shifts an existing atom to a new position between its neighboring atoms.
// If the move crosses matrix elements, accept/reject is based on the resulting
// change in likelihood.
template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::move()
{
    // Select an atom and its left/right neighbors to preserve atomic ordering.
    AtomNeighborhoodType hood = mDomain.randomAtomWithNeighbors(&mRng);
    AtomType *atom = hood.center;
    // Bound the move by neighboring atom positions; use domain endpoints for
    // edge atoms.
    uint64_t lbound = hood.hasLeft() ? hood.left->pos() : 0;
    uint64_t rbound = hood.hasRight() ? hood.right->pos() :
        static_cast<uint64_t>(mDomainLength);

    // Select the new atomic position and map old/new positions to matrix indices.
    uint64_t pos = mRng.uniform64(lbound + 1, rbound - 1);
    unsigned r1 = (atom->pos() / mBinLength) / mNumPatterns;
    unsigned c1 = (atom->pos() / mBinLength) % mNumPatterns;
    unsigned r2 = (pos / mBinLength) / mNumPatterns;
    unsigned c2 = (pos / mBinLength) % mNumPatterns;

    // automatically accept move if it keeps atom in the same matrix element
    if (r1 == r2 && c1 == c2)
    {
        mDomain.move(atom, pos);
        return;
    }
    
    // Conditionally accept moves that change matrix elements.
    float deltaLL = DataModel::deltaLogLikelihood(r1, c1, r2, c2, atom->mass());
    if (std::log(mRng.uniform()) < deltaLL)
    {
        mDomain.move(atom, pos);
        DataModel::safelyChangeMatrix(r1, c1, -atom->mass());
        DataModel::changeMatrix(r2, c2, atom->mass());
    }
}

// Exchange transfers mass between neighboring atoms. When the two atoms affect
// different matrix elements and the conditional Gibbs distribution is defined,
// that distribution proposes the amount of mass to exchange; otherwise the
// exchange is skipped because its conditional variance would be degenerate.
template <class DataModel>
void SingleThreadedGibbsSampler<DataModel>::exchange()
{
    // Select an atom and its right neighbor, wrapping at the end of the domain.
    AtomNeighborhoodType hood = mDomain.randomAtomWithNeighbors(&mRng);
    AtomType *atom1 = hood.center;
    // Exchange with the right neighbor, wrapping the last atom to the first.
    AtomType *atom2 = hood.hasRight() ? hood.right : mDomain.front();

    // Map both atom positions to matrix indices.
    unsigned r1 = (atom1->pos() / mBinLength) / mNumPatterns;
    unsigned c1 = (atom1->pos() / mBinLength) % mNumPatterns;
    unsigned r2 = (atom2->pos() / mBinLength) / mNumPatterns;
    unsigned c2 = (atom2->pos() / mBinLength) % mNumPatterns;

    // Same-bin exchanges only redistribute atomic mass and do not update matrix entries.
    if ((r1 != r2 || c1 != c2) && DataModel::canUseGibbs(c1, c2))
    {
        OptionalFloat mass = DataModel::sampleExchange(r1, c1, atom1->mass(),
            r2, c2, atom2->mass(), &mRng);
        float newMass1 = atom1->mass() + mass.value();
        float newMass2 = atom2->mass() - mass.value();
        if (mass.hasValue() && newMass1 > gaps::epsilon && newMass2 > gaps::epsilon)
        {
            DataModel::safelyChangeMatrix(r1, c1, newMass1 - atom1->mass());
            DataModel::safelyChangeMatrix(r2, c2, newMass2 - atom2->mass());
            atom1->updateMass(newMass1);
            atom2->updateMass(newMass2);
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
