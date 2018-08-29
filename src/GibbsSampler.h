#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#include "AtomicDomain.h"
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
    
    void setSparsity(float alpha, bool singleCell);
    void setMaxGibbsMass(float max);
    void setAnnealingTemp(float temp);
    void setMatrix(const Matrix &mat);

    float chi2() const;
    uint64_t nAtoms() const;

    void sync(const GibbsSampler &sampler);
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

    GapsRng mPropRng;
    AtomicDomain mDomain; // data structure providing access to atoms

    float mAlpha;
    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumPatterns;
    uint64_t mNumBins;
    uint64_t mBinSize;
    uint64_t mDomainLength;

    void makeAndProcessProposal();
    float deathProb(uint64_t nAtoms) const;

    void birth();
    void death();
    void move();
    void exchange();

    void acceptExchange(Atom *a1, Atom *a2, float d1, unsigned r1,
        unsigned c1, unsigned r2, unsigned c2);
    bool updateAtomMass(Atom *atom, float delta);

    OptionalFloat gibbsMass(AlphaParameters alpha, GapsRng *rng);
    OptionalFloat gibbsMass(AlphaParameters alpha, float m1, float m2, GapsRng *rng);

    void safelyChangeMatrix(unsigned row, unsigned col, float delta);
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;

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
mLambda(0.f),
mMaxGibbsMass(100.f),
mAnnealingTemp(1.f),
mNumPatterns(mMatrix.nCol()),
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