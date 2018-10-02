#ifndef __COGAPS_GAPS_RUNNER_H__
#define __COGAPS_GAPS_RUNNER_H__

#include "GapsParameters.h"
#include "GapsResult.h"
#include "GapsStatistics.h"
#include "GibbsSampler.h"

// boost time helpers
#include <boost/date_time/posix_time/posix_time.hpp>
namespace bpt = boost::posix_time;
#define bpt_now() bpt::microsec_clock::local_time()

class GapsRunner
{
public:

    template <class DataType>
    GapsRunner(const DataType &data, const GapsParameters &params);

    template <class DataType>
    void setUncertainty(const DataType &unc, const GapsParameters &params);

    GapsResult run();

private:
    
    GibbsSampler *mASampler;
    GibbsSampler *mPSampler;
    GapsStatistics mStatistics;

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
    void updateSampler(unsigned nA, unsigned nP);
    double estimatedPercentComplete() const;
    void displayStatus();
    void createCheckpoint();
};

template <class DataType>
GapsRunner::GapsRunner(const DataType &data, const GapsParameters &params)
    :
mASampler(new DenseGibbsSampler(data, !params.transposeData, params.nPatterns, params.subsetGenes, params.dataIndicesSubset)),
mPSampler(new DenseGibbsSampler(data, params.transposeData, params.nPatterns, params.subsetGenes, params.dataIndicesSubset)),
mStatistics(mPSampler->dataRows(), mPSampler->dataCols(), params.nPatterns),
mCheckpointOutFile(params.checkpointOutFile),
mMaxIterations(params.nIterations),
mMaxThreads(params.mMaxThreads),
mOutputFrequency(params.mOutputFrequency),
mCheckpointInterval(params.mCheckpointInterval),
mNumPatterns(params.nPatterns),
mNumUpdatesA(0),
mNumUpdatesP(0),
mSeed(params.seed),
mPrintMessages(params.printMessages),
mPrintThreadUsage(params.printThreadUsage),
mPhase('C'),
mFixedMatrix(params.whichFixedMatrix)
{
    mASampler->setSparsity(params.alphaA, params.maxGibbsMassA, params.singleCell);
    mPSampler->setSparsity(params.alphaP, params.maxGibbsMassP, params.singleCell);

    switch (mFixedMatrix)
    {
        case 'A' : mASampler->setMatrix(params.fixedMatrix); break;
        case 'P' : mPSampler->setMatrix(params.fixedMatrix); break;
        default: break;
    }

    // overwrite with info from checkpoint file
    if (params.useCheckPoint)
    {
        Archive ar(params.checkpointFile, ARCHIVE_READ);
        ar >> mNumPatterns >> mSeed >> mMaxIterations >> mFixedMatrix >> mPhase
            >> mCurrentIteration >> mNumUpdatesA >> mNumUpdatesP >> mRng
            >> *mASampler >> *mPSampler;
        GapsRng::load(ar);
    }

    mASampler->sync(mPSampler);
    mPSampler->sync(mASampler);
    mASampler->recalculateAPMatrix();
    mPSampler->recalculateAPMatrix();
}

template <class DataType>
void GapsRunner::setUncertainty(const DataType &unc, const GapsParameters &params)
{
    mASampler->setUncertainty(unc, !params.transposeData, params.nPatterns,
        params.subsetGenes, params.dataIndicesSubset);
    mPSampler->setUncertainty(unc, params.transposeData, params.nPatterns,
        params.subsetGenes, params.dataIndicesSubset);
}

#endif // __COGAPS_GAPS_RUNNER_H__