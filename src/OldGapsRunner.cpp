#include "GapsRunner.h"

#ifdef __GAPS_R_BUILD__
#include <Rcpp.h>
#endif

#ifdef __GAPS_OPENMP__
#include <omp.h>
#endif

///////////////////////////// RAII wrapper /////////////////////////////////////

GapsRunner::~GapsRunner()
{
    delete mRunner;
}

GapsResult GapsRunner::run()
{
    return mRunner->run();
}

///////////////////////// Abstract Interface ///////////////////////////////////

AbstractGapsRunner::AbstractGapsRunner(const GapsParameters &params,
GapsRandomState *randState)
    :
mStatistics(params.nGenes, params.nSamples, params.nPatterns),
mRandState(randState),
mRng(randState),
mCheckpointOutFile(params.checkpointOutFile),
mCurrentIteration(0),
mMaxIterations(params.nIterations),
mMaxThreads(params.maxThreads),
mOutputFrequency(params.outputFrequency),
mCheckpointInterval(params.checkpointInterval),
mNumPatterns(params.nPatterns),
mNumUpdatesA(0),
mNumUpdatesP(0),
mSeed(params.seed),
mPrintMessages(params.printMessages),
mPrintThreadUsage(params.printThreadUsage),
mPhase('C'),
mFixedMatrix(params.whichFixedMatrix)
{}

///////////////////// DenseGapsRunner Implementation ///////////////////////////

void DenseGapsRunner::updateSampler(unsigned nA, unsigned nP)
{
    if (mFixedMatrix != 'A')
    {
        mNumUpdatesA += nA;
        mASampler.update(nA, mMaxThreads);
        if (mFixedMatrix != 'P')
        {
            mPSampler.sync(mASampler, mMaxThreads);
        }
    }

    if (mFixedMatrix != 'P')
    {
        mNumUpdatesP += nP;
        mPSampler.update(nP, mMaxThreads);
        if (mFixedMatrix != 'A')
        {
            mASampler.sync(mPSampler, mMaxThreads);
        }
    }
}

void DenseGapsRunner::setUncertainty(const Matrix &unc, const GapsParameters &params)
{
    mASampler.setUncertainty(unc, !params.transposeData, !params.subsetGenes, params);
    mPSampler.setUncertainty(unc, params.transposeData, params.subsetGenes, params);
}

void DenseGapsRunner::setUncertainty(const std::string &unc, const GapsParameters &params)
{
    mASampler.setUncertainty(unc, !params.transposeData, !params.subsetGenes, params);
    mPSampler.setUncertainty(unc, params.transposeData, params.subsetGenes, params);
}

///////////////////// SparseGapsRunner Implementation //////////////////////////

float SparseGapsRunner::chiSq() const
{
    // doesn't matter which sampler is called
    return mPSampler.chiSq();
}

float SparseGapsRunner::meanChiSq() const
{
    // need to pass P sampler (due to configuration of internal data)
    return mStatistics.meanChiSq(mPSampler);
}

unsigned SparseGapsRunner::nAtoms(char which) const
{
    return which == 'A' ? mASampler.nAtoms() : mPSampler.nAtoms();
}

void SparseGapsRunner::setAnnealingTemp(float temp)
{
    mASampler.setAnnealingTemp(temp);
    mPSampler.setAnnealingTemp(temp);
}

void SparseGapsRunner::updateStatistics()
{
    mStatistics.update(mASampler, mPSampler);
}

Archive& SparseGapsRunner::readSamplers(Archive &ar)
{
    ar >> mASampler >> mPSampler;
    return ar;
}

Archive& SparseGapsRunner::writeSamplers(Archive &ar)
{
    ar << mASampler << mPSampler;
    return ar;
}

void SparseGapsRunner::updateSampler(unsigned nA, unsigned nP)
{
    if (mFixedMatrix != 'A')
    {
        mNumUpdatesA += nA;
        mASampler.update(nA, mMaxThreads);
        if (mFixedMatrix != 'P')
        {
            mPSampler.sync(mASampler, mMaxThreads);
        }
    }

    if (mFixedMatrix != 'P')
    {
        mNumUpdatesP += nP;
        mPSampler.update(nP, mMaxThreads);
        if (mFixedMatrix != 'A')
        {
            mASampler.sync(mPSampler, mMaxThreads);
        }
    }
}

void SparseGapsRunner::setUncertainty(const Matrix &unc, const GapsParameters &params)
{
    // nothing happens - SparseGibbsSampler assumes default uncertainty always
}

void SparseGapsRunner::setUncertainty(const std::string &unc, const GapsParameters &params)
{
    // nothing happens - SparseGibbsSampler assumes default uncertainty always
}
