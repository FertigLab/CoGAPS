#ifndef __COGAPS_ASYNCHRONOUS_GIBBS_SAMPLER_H__
#define __COGAPS_ASYNCHRONOUS_GIBBS_SAMPLER_H__

#include "../atomic/AtomicDomain.h"
#include "../atomic/ProposalQueue.h"
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

//////////////////////////// AsynchronousGibbsSampler Interface ////////////////////////////

class GapsStatistics;

template <class DataModel>
class AsynchronousGibbsSampler;

template <class DataModel>
Archive& operator<<(Archive &ar, const AsynchronousGibbsSampler<DataModel> &sampler);

template <class DataModel>
Archive& operator>>(Archive &ar, AsynchronousGibbsSampler<DataModel> &sampler);

template <class DataModel>
class AsynchronousGibbsSampler : public DataModel
{
public:

    template <class DataType>
    AsynchronousGibbsSampler(const DataType &data, bool transpose, bool subsetRows,
        float alpha, float maxGibbsMass, const GapsParameters &params,
        GapsRandomState *randState);

    unsigned nAtoms() const;
    float getAverageQueueLength() const;

    void update(unsigned nSteps, unsigned nThreads);

    friend Archive& operator<< <DataModel> (Archive &ar, const AsynchronousGibbsSampler &s);
    friend Archive& operator>> <DataModel> (Archive &ar, AsynchronousGibbsSampler &s);

private:

    AtomicDomain mDomain; // data structure providing access to atoms
    ProposalQueue mQueue; // creates queue of proposals that get evaluated by sampler

    float mAvgQueueLength;
    float mNumQueueSamples;

    void birth(const AtomicProposal &prop);
    void death(const AtomicProposal &prop);
    void move(const AtomicProposal &prop);
    void exchange(const AtomicProposal &prop);
};

//////////////////// AsynchronousGibbsSampler - templated functions ////////////////////////

template <class DataModel>
template <class DataType>
AsynchronousGibbsSampler<DataModel>::AsynchronousGibbsSampler(const DataType &data,
bool transpose, bool subsetRows, float alpha, float maxGibbsMass,
const GapsParameters &params, GapsRandomState *randState)
    :
DataModel(data, transpose, subsetRows, params, alpha),
mDomain(DataModel::nElements()),
mQueue(DataModel::nElements(), DataModel::nPatterns(), randState),
mAvgQueueLength(0),
mNumQueueSamples(0)
{
    mQueue.setAlpha(alpha);
    mQueue.setLambda(DataModel::lambda());
}

template <class DataModel>
unsigned AsynchronousGibbsSampler<DataModel>::nAtoms() const
{
    return mDomain.size();
}

template <class DataModel>
float AsynchronousGibbsSampler<DataModel>::getAverageQueueLength() const
{
    return mAvgQueueLength;
}

template <class DataModel>
void AsynchronousGibbsSampler<DataModel>::update(unsigned nSteps, unsigned nThreads)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        // create the largest queue possible, without hitting any conflicts
        mQueue.populate(mDomain, nSteps - n);
        n += mQueue.nProcessed();
    
        if (n < nSteps) // don't count last one since it might be truncated
        {
            mNumQueueSamples += 1.f;
            mAvgQueueLength *= (mNumQueueSamples - 1.f) / mNumQueueSamples;
            mAvgQueueLength += static_cast<float>(mQueue.size()) / mNumQueueSamples;
        }
        
        // process all proposed updates in parallel - the way the queue is 
        // populated ensures no race conditions will happen
        #pragma omp parallel for num_threads(nThreads)
        for (unsigned i = 0; i < mQueue.size(); ++i)
        {
            switch (mQueue[i].type)
            {
                case 'B': birth(mQueue[i]);    break;
                case 'D': death(mQueue[i]);    break;
                case 'M': move(mQueue[i]);     break;
                case 'E': exchange(mQueue[i]); break;
            }
        }
        mQueue.clear();
    }

    GAPS_ASSERT(mDomain.isSorted());
}

// add an atom at a random position, calculate mass either with an
// exponential distribution or with the gibbs mass distribution
template <class DataModel>
void AsynchronousGibbsSampler<DataModel>::birth(const AtomicProposal &prop)
{
    // try to get mass using gibbs, resort to exponential if needed
    OptionalFloat mass = DataModel::canUseGibbs(prop.c1)
        ? DataModel::sampleBirth(prop.r1, prop.c1, &(prop.rng))
        : prop.rng.exponential(DataModel::lambda());

    // accept mass as long as gibbs succeded and it's non-zero
    if (mass.hasValue() && mass.value() >= gaps::epsilon)
    {
        mQueue.acceptBirth();
        prop.atom1->mass = mass.value();
        DataModel::changeMatrix(prop.r1, prop.c1, mass.value());
        return;
    }

    // otherwise reject birth
    mQueue.rejectBirth();
    mDomain.erase(prop.atom1->pos);
}

// attempt to rebirth an atom in place of the killed atom
template <class DataModel>
void AsynchronousGibbsSampler<DataModel>::death(const AtomicProposal &prop)
{
    // try to do a rebirth in the place of this atom
    if (DataModel::canUseGibbs(prop.c1))
    {
        OptionalFloat mass = DataModel::sampleDeathAndRebirth(prop.r1, prop.c1,
            -1.f * prop.atom1->mass, &(prop.rng));
        if (mass.hasValue() && mass.value() >= gaps::epsilon)
        {
            mQueue.rejectDeath();
            DataModel::safelyChangeMatrix(prop.r1, prop.c1, mass.value() - prop.atom1->mass);
            prop.atom1->mass = mass.value();
            return;
        }
    }

    // if rebirth fails, then kill off atom
    mQueue.acceptDeath();
    DataModel::safelyChangeMatrix(prop.r1, prop.c1, -1.f * prop.atom1->mass);
    mDomain.erase(prop.atom1->pos);
}

// move mass from src to dest in the atomic domain
template <class DataModel>
void AsynchronousGibbsSampler<DataModel>::move(const AtomicProposal &prop)
{
    GAPS_ASSERT(prop.r1 != prop.r2 || prop.c1 != prop.c2);

    float deltaLL = DataModel::deltaLogLikelihood(prop.r1, prop.c1, prop.r2, prop.c2,
        prop.atom1->mass);
    if (std::log(prop.rng.uniform()) < deltaLL)
    {
        prop.atom1->pos = prop.pos;
        DataModel::safelyChangeMatrix(prop.r1, prop.c1, -prop.atom1->mass);
        DataModel::changeMatrix(prop.r2, prop.c2, prop.atom1->mass);
        return;
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small
template <class DataModel>
void AsynchronousGibbsSampler<DataModel>::exchange(const AtomicProposal &prop)
{
    GAPS_ASSERT(prop.r1 != prop.r2 || prop.c1 != prop.c2);

    if (DataModel::canUseGibbs(prop.c1, prop.c2))
    {
        OptionalFloat mass = DataModel::sampleExchange(prop.r1, prop.c1, prop.atom1->mass,
            prop.r2, prop.c2, prop.atom2->mass, &(prop.rng));
        float newMass1 = prop.atom1->mass + mass.value();
        float newMass2 = prop.atom2->mass - mass.value();
        if (mass.hasValue() && newMass1 > gaps::epsilon && newMass2 > gaps::epsilon)
        {
            prop.atom1->mass = newMass1;
            prop.atom2->mass = newMass2;
            DataModel::safelyChangeMatrix(prop.r1, prop.c1, newMass1 - prop.atom1->mass);
            DataModel::safelyChangeMatrix(prop.r2, prop.c2, newMass2 - prop.atom2->mass);
            return;
        }
    }
}

template <class DataModel>
Archive& operator<<(Archive &ar, const AsynchronousGibbsSampler<DataModel> &s)
{
    operator<<(ar, static_cast<const DataModel&>(s)) << s.mDomain << s.mQueue;
    return ar;
}

template <class DataModel>
Archive& operator>>(Archive &ar, AsynchronousGibbsSampler<DataModel> &s)
{
    operator>>(ar, static_cast<DataModel&>(s)) >> s.mDomain >> s.mQueue;
    return ar;
}

#endif // __COGAPS_ASYNCHRONOUS_GIBBS_SAMPLER_H__