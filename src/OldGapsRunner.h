#ifndef __COGAPS_GAPS_RUNNER_H__
#define __COGAPS_GAPS_RUNNER_H__

#include "GapsParameters.h"
#include "GapsResult.h"
#include "GapsStatistics.h"
#include "gibbs_sampler/GibbsSampler.h"
#include "gibbs_sampler/DenseGibbsSampler.h"
#include "gibbs_sampler/SparseGibbsSampler.h"

#include <string>

// boost time helpers
#include <boost/date_time/posix_time/posix_time.hpp>
namespace bpt = boost::posix_time;
#define bpt_now() bpt::microsec_clock::local_time()

// forward declarations
class AbstractGapsRunner;

///////////////////////////// RAII wrapper /////////////////////////////////////

// This is the class that is exposed to the top-level CoGAPS routine - all 
// aspects of CoGAPS can be managed through this class. The class itself is
// just a lightweight wrapper around an abstract interface, which allows for
// multiple types of GapsRunner to be declared. Which implementation is used
// depends on the parameters passed to the GapsRunner constructor.
class GapsRunner
{
public:

    template <class DataType>
    GapsRunner(const DataType &data, const GapsParameters &params,
        GapsRandomState *randState);

    ~GapsRunner();

    template <class DataType>
    void setUncertainty(const DataType &unc, const GapsParameters &params);

    GapsResult run();

private:

    AbstractGapsRunner *mRunner;

    GapsRunner(const GapsRunner &p); // don't allow copies
    GapsRunner& operator=(const GapsRunner &p); // don't allow copies    
};

///////////////////////// Abstract Interface ///////////////////////////////////

// This class is the abstract interface that any implementation of GapsRunner
// must satisfy. It provides a factory method that will create the appropiate
// derived class depending on the parameters passed in.
class AbstractGapsRunner
{
public:

    AbstractGapsRunner(const GapsParameters &params, GapsRandomState *randState);
    virtual ~AbstractGapsRunner() {}

    template <class DataType>
    static AbstractGapsRunner* create(const DataType &data,
        const GapsParameters &params, GapsRandomState *randState);

    // can't use template with virtual function
    virtual void setUncertainty(const Matrix &unc, const GapsParameters &params) = 0;
    virtual void setUncertainty(const std::string &unc, const GapsParameters &params) = 0;

    GapsResult run();

protected:

    GapsStatistics mStatistics;

    const GapsRandomState *mRandState; // used for writing state to checkpoint
    mutable GapsRng mRng;

    std::string mCheckpointOutFile;

    bpt::ptime mStartTime;

    unsigned mCurrentIteration;
    unsigned mMaxIterations;
    unsigned mMaxThreads;
    unsigned mOutputFrequency;
    unsigned mCheckpointInterval;
    unsigned mNumPatterns;
    unsigned mNumUpdatesA;
    unsigned mNumUpdatesP;
    uint32_t mSeed;

    bool mPrintMessages;
    bool mPrintThreadUsage;

    char mPhase;
    char mFixedMatrix;
        
    void runOnePhase();
    double estimatedPercentComplete() const;
    void displayStatus();
    void createCheckpoint();

    virtual float chiSq() const = 0;
    virtual float meanChiSq() const = 0;
    virtual unsigned nAtoms(char which) const = 0;
    virtual void setAnnealingTemp(float temp) = 0;
    virtual void updateStatistics() = 0;
    virtual Archive& readSamplers(Archive &ar) = 0;
    virtual Archive& writeSamplers(Archive &ar) = 0;
    virtual void updateSampler(unsigned nA, unsigned nP) = 0;
};

///////////////////// GapsRunner Implementations ///////////////////////////////

// This implementation uses a DenseGibbsSampler internally
class DenseGapsRunner : public AbstractGapsRunner
{
public:

    ~DenseGapsRunner() {}

    template <class DataType>
    DenseGapsRunner(const DataType &data, const GapsParameters &params,
        GapsRandomState *randState);

    void setUncertainty(const Matrix &unc, const GapsParameters &params);
    void setUncertainty(const std::string &unc, const GapsParameters &params);

private:

    DenseGibbsSampler mASampler;
    DenseGibbsSampler mPSampler;

    float chiSq() const;
    float meanChiSq() const;
    unsigned nAtoms(char which) const;
    void setAnnealingTemp(float temp);
    void updateStatistics();
    Archive& readSamplers(Archive &ar);
    Archive& writeSamplers(Archive &ar);
    void updateSampler(unsigned nA, unsigned nP);
};

// This implementation uses a SparseGibbsSampler internally
class SparseGapsRunner : public AbstractGapsRunner
{
public:

    ~SparseGapsRunner() {}

    template <class DataType>
    SparseGapsRunner(const DataType &data, const GapsParameters &params,
        GapsRandomState *randState);

    void setUncertainty(const Matrix &unc, const GapsParameters &params);
    void setUncertainty(const std::string &unc, const GapsParameters &params);

private:

    SparseGibbsSampler mASampler;
    SparseGibbsSampler mPSampler;

    float chiSq() const;
    float meanChiSq() const;
    unsigned nAtoms(char which) const;
    void setAnnealingTemp(float temp);
    void updateStatistics();
    Archive& readSamplers(Archive &ar);
    Archive& writeSamplers(Archive &ar);
    void updateSampler(unsigned nA, unsigned nP);
};

/////////////////////// GapsRunner - templated functions ///////////////////////

template <class DataType>
GapsRunner::GapsRunner(const DataType &data, const GapsParameters &params,
GapsRandomState *randState)
    : mRunner(AbstractGapsRunner::create(data, params, randState))
{}

template <class DataType>
void GapsRunner::setUncertainty(const DataType &unc, const GapsParameters &params)
{
    mRunner->setUncertainty(unc, params);
}

/////////////////// AbstractGapsRunner - templated functions ///////////////////

template <class DataType>
AbstractGapsRunner* AbstractGapsRunner::create(const DataType &data,
const GapsParameters &params, GapsRandomState *randState)
{
    if (params.useSparseOptimization)
    {
        return new SparseGapsRunner(data, params, randState);
    }
    return new DenseGapsRunner(data, params, randState);
}

//////////////////// DenseGapsRunner - templated functions /////////////////////

template <class DataType>
DenseGapsRunner::DenseGapsRunner(const DataType &data,
const GapsParameters &params, GapsRandomState *randState)
    :
AbstractGapsRunner(params, randState),
mASampler(data, !params.transposeData, !params.subsetGenes, params.alphaA,
    params.maxGibbsMassA, params, randState),
mPSampler(data, params.transposeData, params.subsetGenes, params.alphaP,
    params.maxGibbsMassP, params, randState)
{
    switch (mFixedMatrix)
    {
        case 'A' : mASampler.setMatrix(params.fixedMatrix); break;
        case 'P' : mPSampler.setMatrix(params.fixedMatrix); break;
        default: break; // 'N' for none
    }

    // overwrite with info from checkpoint file
    if (params.useCheckPoint)
    {
        Archive ar(params.checkpointFile, ARCHIVE_READ);
        ar >> mNumPatterns >> mSeed >> mMaxIterations >> mFixedMatrix >> mPhase
            >> mCurrentIteration >> mNumUpdatesA >> mNumUpdatesP >> mRng;
        readSamplers(ar);
        randState->load(ar);
    }

    mPSampler.sync(mASampler);
    mASampler.sync(mPSampler);

    // AP matrix not stored in checkpoint
    if (params.useCheckPoint)
    {
        mASampler.recalculateAPMatrix();
        mPSampler.recalculateAPMatrix();
    }
}

//////////////////// SparseGapsRunner - templated functions ////////////////////

template <class DataType>
SparseGapsRunner::SparseGapsRunner(const DataType &data,
const GapsParameters &params, GapsRandomState *randState)
    :
AbstractGapsRunner(params, randState),
mASampler(data, !params.transposeData, !params.subsetGenes, params.alphaA,
    params.maxGibbsMassA, params, randState),
mPSampler(data, params.transposeData, params.subsetGenes, params.alphaP,
    params.maxGibbsMassP, params, randState)
{
    switch (mFixedMatrix)
    {
        case 'A' : mASampler.setMatrix(params.fixedMatrix); break;
        case 'P' : mPSampler.setMatrix(params.fixedMatrix); break;
        default: break;
    }

    // overwrite with info from checkpoint file
    if (params.useCheckPoint)
    {
        Archive ar(params.checkpointFile, ARCHIVE_READ);
        ar >> mNumPatterns >> mSeed >> mMaxIterations >> mFixedMatrix >> mPhase
            >> mCurrentIteration >> mNumUpdatesA >> mNumUpdatesP >> mRng;
        readSamplers(ar);
        randState->load(ar);
    }

    mASampler.sync(mPSampler);
    mPSampler.sync(mASampler);
}

#endif // __COGAPS_GAPS_RUNNER_H__
