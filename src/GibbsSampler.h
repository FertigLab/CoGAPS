#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#include "AtomicDomain.h"
#include "ProposalQueue.h"
#include "data_structures/Matrix.h"
#include "math/Algorithms.h"
#include "math/Random.h"

#include <vector>

class GibbsSampler
{
public:

    template <class DataType>
    GibbsSampler(const DataType &data, bool transpose, unsigned nPatterns,
        bool partitionRows, const std::vector<unsigned> &indices);

    template <class DataType>
    void setUncertainty(const DataType &unc, bool transposeData,
        bool partitionRows, const std::vector<unsigned> &indices);

    unsigned dataRows() const;
    unsigned dataCols() const;
    
    void setSparsity(float alpha, float maxGibbsMass, bool singleCell);
    void setAnnealingTemp(float temp);
    void setMatrix(const Matrix &mat);

    float chi2() const;
    uint64_t nAtoms() const;

    void recalculateAPMatrix();
    void sync(const GibbsSampler &sampler, unsigned nThreads=1);
    void update(unsigned nSteps, unsigned nCores);

    // serialization
    friend Archive& operator<<(Archive &ar, GibbsSampler &sampler);
    friend Archive& operator>>(Archive &ar, GibbsSampler &sampler);

    #ifdef GAPS_DEBUG
    bool internallyConsistent();
    #endif

private:

    friend class GapsStatistics;

    ColMatrix mDMatrix; // samples by genes for A, genes by samples for P
    ColMatrix mSMatrix; // same configuration as D
    ColMatrix mAPMatrix; // cached product of A and P, same configuration as D

    ColMatrix mMatrix; // genes by patterns for A, samples by patterns for P
    const ColMatrix *mOtherMatrix; // pointer to P if this is A, and vice versa

    AtomicDomain mDomain; // data structure providing access to atoms
    ProposalQueue mQueue; // creates queue of proposals that get evaluated by sampler

    float mAlpha;
    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumPatterns;
    uint64_t mNumBins;
    uint64_t mBinLength;
    uint64_t mDomainLength;

    void processProposal(const AtomicProposal &prop);
    float deathProb(uint64_t nAtoms) const;

    void birth(const AtomicProposal &prop);
    void death(const AtomicProposal &prop);
    void move(const AtomicProposal &prop);
    void exchange(const AtomicProposal &prop);
    void exchangeUsingMetropolisHastings(const AtomicProposal &prop,
        AlphaParameters alpha);
    void acceptExchange(const AtomicProposal &prop, float delta);
    bool updateAtomMass(Atom *atom, float delta);

    void changeMatrix(unsigned row, unsigned col, float delta);
    void safelyChangeMatrix(unsigned row, unsigned col, float delta);
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    bool canUseGibbs(unsigned col) const;
    bool canUseGibbs(unsigned c1, unsigned c2) const;

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    AlphaParameters alphaParametersWithChange(unsigned row, unsigned col, float ch);
};

template <class DataType>
GibbsSampler::GibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool partitionRows,
const std::vector<unsigned> &indices)
    :
mDMatrix(data, transposeData, partitionRows, indices),
mSMatrix(gaps::algo::pmax(mDMatrix, 0.1f)),
mAPMatrix(mDMatrix.nRow(), mDMatrix.nCol()),
mMatrix(mDMatrix.nCol(), nPatterns),
mOtherMatrix(NULL),
mDomain(mMatrix.nRow() * mMatrix.nCol()),
mQueue(mMatrix.nRow(), mMatrix.nCol()),
mLambda(0.f),
mMaxGibbsMass(0.f),
mAnnealingTemp(1.f),
mNumPatterns(mMatrix.nCol()),
mNumBins(mMatrix.nRow() * mMatrix.nCol()),
mBinLength(std::numeric_limits<uint64_t>::max() / mNumBins),
mDomainLength(mBinLength * mNumBins)
{
    // default sparsity parameters
    setSparsity(0.01, 100.f, false);
}

template <class DataType>
void GibbsSampler::setUncertainty(const DataType &unc,
bool transpose, bool partitionRows, const std::vector<unsigned> &indices)
{
    mSMatrix = ColMatrix(unc, transpose, partitionRows, indices);
}

#endif // __COGAPS_GIBBS_SAMPLER_H__