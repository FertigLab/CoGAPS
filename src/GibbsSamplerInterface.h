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

// This is a polymorphic wrapper to an underlying implementation. The purpose of
// this is to move the virtual functions to the highest level possible where
// they are called infrequently, rather than pay the performance cost of
// dynamic dispatch of frequently called, inexpensive functions.

// interface
class GibbsSampler
{
public:

    virtual unsigned dataRows() const = 0;
    virtual unsigned dataCols() const = 0;
    
    virtual void setSparsity(float alpha, float maxGibbsMass, bool singleCell) = 0;
    virtual void setAnnealingTemp(float temp) = 0;
    virtual void setMatrix(const Matrix &mat) = 0;

    virtual float chi2() const = 0;
    virtual uint64_t nAtoms() const = 0;

    virtual void recalculateAPMatrix() = 0;
    virtual void sync(const GibbsSampler *sampler, unsigned nThreads=1) = 0;
    virtual void update(unsigned nSteps, unsigned nCores) = 0;
};

// wrapper for a dense GibbsSampler implementation - all data in this class
// is stored as a dense matrix
class DenseGibbsSampler : public GibbsSampler
{
public:

    template <class DataType>
    DenseGibbsSampler(const DataType &data, bool transpose, unsigned nPatterns,
    bool partitionRows, const std::vector<unsigned> &indices)
        : mImplementation(data, transpose, nPatterns, partitionRows, indices);
    {}

    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, unsigned nPatterns,
    bool partitionRows, const std::vector<unsigned> &indices)
    {
        mImplementation.setUncertainty(data, transpose, nPatterns, partitionRows, indices);
    }

    unsigned dataRows() const { return mImplementation.dataRows(); }
    unsigned dataCols() const { return mImplementation.dataCols(); }

    void setSparsity(float alpha, float maxGibbsMass, bool singleCell) { mImplementation.setSparsity(alpha, maxGibbsMass, singleCell); }
    void setAnnealingTemp(float temp) { mImplementation.setAnnealingTemp(temp); }
    void setMatrix(const Matrix &mat) { mImplementation.setMatrix(mat); }

    float chi2() const { return mImplementation.chi2(); }
    uint64_t nAtoms() const { return mImplementation.nAtoms(); }

    void recalculateAPMatrix() { mImplementation.recalculateAPMatrix(); }
    void sync(const GibbsSampler *sampler, unsigned nThreads=1) { mImplementation.sync(sampler->mImplementation, nThreads); }
    void update(unsigned nSteps, unsigned nCores) { mImplementation.update(nSteps, nCores); }

private:

    DenseGibbsSamplerImplementation mImplementation;
};

// wrapper for a sparse GibbsSampler implementation - all data in this class
// is stored as a sparse matrix
class SparseGibbsSampler : public GibbsSampler
{
public:

    template <class DataType>
    SparseGibbsSampler(const DataType &data, bool transpose, unsigned nPatterns,
    bool partitionRows, const std::vector<unsigned> &indices)
        : mImplementation(data, transpose, nPatterns, partitionRows, indices);
    {}

    template <class DataType>
    void setUncertainty(const DataType &data, bool transpose, unsigned nPatterns,
    bool partitionRows, const std::vector<unsigned> &indices)
    {
        GAPS_ASSERT(false); // should never reach
    }

    unsigned dataRows() const { return mImplementation.dataRows(); }
    unsigned dataCols() const { return mImplementation.dataCols(); }

    void setSparsity(float alpha, float maxGibbsMass, bool singleCell) { mImplementation.setSparsity(alpha, maxGibbsMass, singleCell); }
    void setAnnealingTemp(float temp) { mImplementation.setAnnealingTemp(temp); }
    void setMatrix(const Matrix &mat) { mImplementation.setMatrix(mat); }

    float chi2() const { return mImplementation.chi2(); }
    uint64_t nAtoms() const { return mImplementation.nAtoms(); }

    void recalculateAPMatrix() { GAPS_ASSERT(false); /* should never reach */}
    void sync(const GibbsSampler &sampler, unsigned nThreads=1) { mImplementation.sync(sampler); }
    void update(unsigned nSteps, unsigned nCores) { mImplementation.update(nSteps, nCores); }

private:
    
    SparseGibbsSamplerImplementation mImplementation;
};

#endif // __COGAPS_GIBBS_SAMPLER_H__