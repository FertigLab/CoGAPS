#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#include "AtomicDomain.h"
#include "ProposalQueue.h"
#include "data_structures/Matrix.h"
#include "math/Algorithms.h"
#include "math/Random.h"
#include "math/SIMD.h"
#include "utils/Archive.h"
#include "utils/GapsAssert.h"

#include <algorithm>

class GibbsSampler
{
private:

    friend class GapsStatistics;

    ColMatrix mDMatrix;
    ColMatrix mSMatrix;
    ColMatrix mAPMatrix;

    ColMatrix mMatrix;
    const ColMatrix* mOtherMatrix;

    GapsRng mPropRng;
    AtomicDomain mDomain;

    float mAlpha;
    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumRows;
    unsigned mNumCols;

    uint64_t mNumBins;
    uint64_t mBinSize;
    uint64_t mDomainLength;

    float mAvgQueue;
    float mNumQueues;

    void makeAndProcessProposal();
    float deathProb(uint64_t nAtoms) const;

    void birth();
    void death();
    void move();
    void exchange();

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;

    bool canUseGibbs(unsigned col) const;
    bool canUseGibbs(unsigned c1, unsigned c2) const;

    void updateAPMatrix(unsigned row, unsigned col, float delta);

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2,
        unsigned c2);

    float computeDeltaLL(unsigned row, unsigned col, float mass);
    float computeDeltaLL(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);

    bool updateAtomMass(Atom *atom, float delta);
    void acceptExchange(AtomicProposal *prop, float d1, unsigned r1,
        unsigned c1, unsigned r2, unsigned c2);

    OptionalFloat gibbsMass(AlphaParameters alpha, GapsRng *rng);
    OptionalFloat gibbsMass(AlphaParameters alpha, float m1, float m2, GapsRng *rng);

public:

    template <class DataType>
    GibbsSampler(const DataType &data, bool transpose, unsigned nPatterns,
        bool partitionRows, const std::vector<unsigned> &indices);

    template <class DataType>
    void setUncertainty(const DataType &unc, bool transposeData,
        bool partitionRows, const std::vector<unsigned> &indices);
    
    void setSparsity(float alpha, bool singleCell);
    void setMaxGibbsMass(float max);
    void setAnnealingTemp(float temp);

    void setMatrix(const Matrix &mat);

    void update(unsigned nSteps, unsigned nCores);

    unsigned dataRows() const;
    unsigned dataCols() const;

    float chi2() const;
    uint64_t nAtoms() const;

    void sync(GibbsSampler &sampler);

    #ifdef GAPS_DEBUG
    bool internallyConsistent();
    #endif

    // serialization
    friend Archive& operator<<(Archive &ar, GibbsSampler &sampler);
    friend Archive& operator>>(Archive &ar, GibbsSampler &sampler);
};

template <class DataType>
GibbsSampler::GibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool partitionRows,
const std::vector<unsigned> &indices)
    :
mDMatrix(data, transposeData, partitionRows, indices),
mSMatrix(mDMatrix.pmax(0.1f, 0.1f)),
mAPMatrix(mDMatrix.nRow(), mDMatrix.nCol()),
mMatrix(mDMatrix.nCol(), nPatterns),
mOtherMatrix(NULL),
mDomain(mMatrix.nRow() * mMatrix.nCol()),
mLambda(0.f),
mMaxGibbsMass(100.f),
mAnnealingTemp(1.f),
mNumRows(mMatrix.nRow()),
mNumCols(mMatrix.nCol()),
mAvgQueue(0.f),
mNumQueues(0.f),
mNumBins(mMatrix.nRow() * mMatrix.nCol()),
mBinSize(std::numeric_limits<uint64_t>::max() / mNumBins),
mDomainLength(mBinSize * mNumBins)
{
    // default sparsity parameters
    setSparsity(0.01, false);
}

template <class DataType>
void GibbsSampler::setUncertainty(const DataType &unc,
bool transpose, bool partitionRows, const std::vector<unsigned> &indices)
{
    mSMatrix = ColMatrix(unc, transpose, partitionRows, indices);
}

#endif // __COGAPS_GIBBS_SAMPLER_H__