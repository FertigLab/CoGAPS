#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#define DEFAULT_ALPHA           0.01f
#define DEFAULT_MAX_GIBBS_MASS  100.f

#include "AtomicDomain.h"
#include "ProposalQueue.h"
#include "data_structures/Matrix.h"
#include "math/Algorithms.h"
#include "math/Random.h"

#include <vector>

////////////////////////////// CLASS DEFINITIONS ///////////////////////////////

// These classes provide the various implementations of a GibbsSampler. Compile
// time polymorphism is used to reduce code duplication.

// can't be constructed, interface is avaiable through derived classes
template <class T>
class GibbsSamplerImplementation
{
private:

    friend class T;         

    GibbsSamplerImplementation();

    void setAnnealingTemp(float temp);

    void update(unsigned nSteps, unsigned nThreads);
    void processProposal(const AtomicProposal &prop);
    void birth(const AtomicProposal &prop);
    void death(const AtomicProposal &prop);
    void move(const AtomicProposal &prop);
    void exchange(const AtomicProposal &prop);
    void exchangeUsingMetropolisHastings(const AtomicProposal &prop,
        AlphaParameters alpha);
    void acceptExchange(const AtomicProposal &prop, float delta);
    bool updateAtomMass(Atom *atom, float delta);

    AtomicDomain mDomain; // data structure providing access to atoms
    ProposalQueue mQueue; // creates queue of proposals that get evaluated by sampler

    float mAlpha;
    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumPatterns;
    uint64_t mBinLength;
};

class DenseGibbsSamplerImplementation : public GibbsSamplerImplementation<DenseGibbsSamplerImplementation>
{
private:

    friend class DenseGibbsSampler;

    DenseColMatrix mDMatrix; // samples by genes for A, genes by samples for P
    DenseColMatrix mSMatrix; // same configuration as D
    DenseColMatrix mAPMatrix; // cached product of A and P, same configuration as D
    DenseColMatrix mMatrix; // genes by patterns for A, samples by patterns for P
    const DenseColMatrix *mOtherMatrix; // pointer to P if this is A, and vice versa
    
    // private constructor allows only DenseGibbsSampler to construct this
    DenseGibbsSamplerImplementation(const DataType &data, bool transpose,
        unsigned nPatterns, bool partitionRows,
        const std::vector<unsigned> &indices);

    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, unsigned nPatterns,
        bool partitionRows, const std::vector<unsigned> &indices);

    unsigned dataRows() const;
    unsigned dataCols() const;

    void setSparsity(float alpha, float maxGibbsMass, bool singleCell);
    void setMatrix(const Matrix &mat);

    void sync(const GibbsSampler &sampler)

    void changeMatrix(unsigned row, unsigned col, float delta);
    void safelyChangeMatrix(unsigned row, unsigned col, float delta);
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    bool canUseGibbs(unsigned col) const;
    bool canUseGibbs(unsigned c1, unsigned c2) const;
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);
};

class SparseGibbsSamplerImplementation : public GibbsSamplerImplementation<SparseGibbsSamplerImplementation>
{
private :

    friend class SparseGibbsSampler;

    SparseColMatrix mDMatrix;
    HybridColMatrix mSparseMatrix;
    DenseRowMatrix mDenseMatrix;
    const HybridColMatrix *mOtherSparseMatrix;
    const DenseRowMatrix *mOtherDenseMatrix;

    // private constructor allows only SparseGibbsSampler to construct this
    SparseGibbsSamplerImplementation(const DataType &data, bool transpose,
        unsigned nPatterns, bool partitionRows,
        const std::vector<unsigned> &indices);

    unsigned dataRows() const;
    unsigned dataCols() const;

    void setSparsity(float alpha, float maxGibbsMass, bool singleCell);
    void setMatrix(const Matrix &mat);

    float chiSq() const;

    bool canUseGibbs(unsigned col) const;
    bool canUseGibbs(unsigned c1, unsigned c2) const;
    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);
};

//////////////////////IMPLEMENTATION OF TEMPLATED FUNCTIONS ////////////////////

template <class T>
void GibbsSamplerImplementation<T>::update(unsigned nSteps, unsigned nThreads)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        // populate queue, prepare domain for this queue
        mQueue.populate(mDomain, nSteps - n);
        n += mQueue.size();
        
        // process all proposed updates
        #pragma omp parallel for num_threads(nThreads)
        for (unsigned i = 0; i < mQueue.size(); ++i)
        {
            processProposal(mQueue[i]);
        }
        mQueue.clear();
    }

    GAPS_ASSERT(internallyConsistent());
    GAPS_ASSERT(mDomain.isSorted());
}

template <class T>
void GibbsSamplerImplementation<T>::processProposal(const AtomicProposal &prop)
{
    switch (prop.type)
    {
        case 'B':
            birth(prop);
            break;
        case 'D':
            death(prop);
            break;
        case 'M':
            move(prop);
            break;
        case 'E':
            exchange(prop);
            break;
    }
}

// add an atom at a random position, calculate mass either with an
// exponential distribution or with the gibbs mass distribution
void GibbsSampler::birth(const AtomicProposal &prop)
{
    // calculate proposed mass
    float mass = canUseGibbs(prop.c1)
        ? gibbsMass(alphaParameters(prop.r1, prop.c1), &(prop.rng)).value()
        : prop.rng.exponential(mLambda);

    // accept mass as long as it's non-zero
    if (mass >= gaps::epsilon)
    {
        mQueue.acceptBirth();
        prop.atom1->mass = mass;
        changeMatrix(prop.r1, prop.c1, mass);
    }
    else
    {
        mQueue.rejectBirth();
        mDomain.erase(prop.atom1->pos);
    }
}






#endif // __COGAPS_GIBBS_SAMPLER_H__