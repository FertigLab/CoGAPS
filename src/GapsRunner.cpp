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

AbstractGapsRunner::AbstractGapsRunner(const GapsParameters &params)
    :
mStatistics(params.nGenes, params.nSamples, params.nPatterns),
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

GapsResult AbstractGapsRunner::run()
{
    GAPS_ASSERT(mPhase == 'C' || mPhase == 'S');

    mStartTime = bpt_now();

    // check if running in debug mode
    #ifdef GAPS_DEBUG
    gaps_printf("Running in debug mode\n");
    #endif

    // calculate appropiate number of threads if compiled with openmp
    #ifdef __GAPS_OPENMP__
    if (mPrintMessages && mPrintThreadUsage)
    {
        unsigned availableThreads = omp_get_max_threads();
        mMaxThreads = gaps::min(availableThreads, mMaxThreads);
        gaps_printf("Running on %d out of %d available threads\n",
            mMaxThreads, availableThreads);
    }
    #endif

    // cascade through phases, allows algorithm to be resumed in either phase
    switch (mPhase)
    {
        case 'C':
            if (mPrintMessages)
            {
                gaps_printf("-- Calibration Phase --\n");
            }
            runOnePhase();
            mPhase = 'S';
            mCurrentIteration = 0;

        case 'S':
            if (mPrintMessages)
            {
                gaps_printf("-- Sampling Phase --\n");
            }
            runOnePhase();
            break;
    }
    GapsResult result(mStatistics);
    result.meanChiSq = meanChiSq();
    return result;    
}

void AbstractGapsRunner::runOnePhase()
{
    for (; mCurrentIteration < mMaxIterations; ++mCurrentIteration)
    {
        createCheckpoint();

        #ifdef __GAPS_R_BUILD__
        Rcpp::checkUserInterrupt();
        #endif
        
        // set annealing temperature in calibration phase
        if (mPhase == 'C')
        {        
            float temp = static_cast<float>(2 * mCurrentIteration)
                / static_cast<float>(mMaxIterations);
            setAnnealingTemp(gaps::min(1.f, temp));
        }
    
        // number of updates per iteration is poisson 
        unsigned nA = mRng.poisson(gaps::max(nAtoms('A'), 10));
        unsigned nP = mRng.poisson(gaps::max(nAtoms('P'), 10));
        updateSampler(nA, nP);

        if (mPhase == 'S')
        {
            updateStatistics();
        }
        displayStatus();
    }
}

// sum coef * log(i) for i = 1 to total, fit coef from number of atoms
// approximates sum of number of atoms (stirling approx to factorial)
// this should be proportional to total number of updates
static double estimatedNumUpdates(double current, double total, float nAtoms)
{
    double coef = nAtoms / std::log(current);
    return coef * std::log(std::sqrt(2 * total * gaps::pi)) +
        total * coef * std::log(total) - total * coef;
}


double AbstractGapsRunner::estimatedPercentComplete() const
{
    double nIter = static_cast<double>(mCurrentIteration);
    double nAtomsA = static_cast<double>(nAtoms('A'));
    double nAtomsP = static_cast<double>(nAtoms('P'));
    
    if (mPhase == 'S')
    {
        nIter += mMaxIterations;
    }

    double totalIter = 2.0 * static_cast<double>(mMaxIterations);

    double estimatedCompleted = estimatedNumUpdates(nIter, nIter, nAtomsA) + 
        estimatedNumUpdates(nIter, nIter, nAtomsP);

    double estimatedTotal = estimatedNumUpdates(nIter, totalIter, nAtomsA) + 
        estimatedNumUpdates(nIter, totalIter, nAtomsP);

    return estimatedCompleted / estimatedTotal;
}

void AbstractGapsRunner::displayStatus()
{
    if (mPrintMessages && mOutputFrequency > 0 && ((mCurrentIteration + 1) % mOutputFrequency) == 0)
    {
        bpt::time_duration diff = bpt_now() - mStartTime;
        double nSecondsCurrent = diff.total_seconds();
        double nSecondsTotal = nSecondsCurrent / estimatedPercentComplete();

        unsigned elapsedSeconds = static_cast<unsigned>(nSecondsCurrent);
        unsigned totalSeconds = static_cast<unsigned>(nSecondsTotal);

        unsigned elapsedHours = elapsedSeconds / 3600;
        elapsedSeconds -= elapsedHours * 3600;
        unsigned elapsedMinutes = elapsedSeconds / 60;
        elapsedSeconds -= elapsedMinutes * 60;

        unsigned totalHours = totalSeconds / 3600;
        totalSeconds -= totalHours * 3600;
        unsigned totalMinutes = totalSeconds / 60;
        totalSeconds -= totalMinutes * 60;

        gaps_printf("%d of %d, Atoms: %lu(%lu), ChiSq: %.0f, Time: %02d:%02d:%02d / %02d:%02d:%02d\n",
            mCurrentIteration + 1, mMaxIterations, nAtoms('A'),
            nAtoms('P'), chiSq(), elapsedHours, elapsedMinutes,
            elapsedSeconds, totalHours, totalMinutes, totalSeconds);
        gaps_flush();
    }
}

void AbstractGapsRunner::createCheckpoint()
{
    if (mCheckpointInterval > 0 && ((mCurrentIteration + 1) % mCheckpointInterval) == 0)
    {
        // create backup file
        std::rename(mCheckpointOutFile.c_str(), (mCheckpointOutFile + ".backup").c_str());
    
        // create checkpoint file
        Archive ar(mCheckpointOutFile, ARCHIVE_WRITE);
        ar << mNumPatterns << mSeed << mMaxIterations << mFixedMatrix << mPhase
            << mCurrentIteration << mNumUpdatesA << mNumUpdatesP << mRng;
        writeSamplers(ar);
        GapsRng::save(ar);

        // delete backup file
        std::remove((mCheckpointOutFile + ".backup").c_str());
    }
}

///////////////////// DenseGapsRunner Implementation ///////////////////////////

float DenseGapsRunner::chiSq() const
{
    // doesn't matter which sampler is called
    return mPSampler.chiSq();
}

float DenseGapsRunner::meanChiSq() const
{
    // need to pass P sampler (due to configuration of internal data)
    return mStatistics.meanChiSq(mPSampler);
}

unsigned DenseGapsRunner::nAtoms(char which) const
{
    return which == 'A' ? mASampler.nAtoms() : mPSampler.nAtoms();
}

void DenseGapsRunner::setAnnealingTemp(float temp)
{
    mASampler.setAnnealingTemp(temp);
    mPSampler.setAnnealingTemp(temp);
}

void DenseGapsRunner::updateStatistics()
{
    mStatistics.update(mASampler, mPSampler);
}

Archive& DenseGapsRunner::readSamplers(Archive &ar)
{
    ar >> mASampler >> mPSampler;
    return ar;
}

Archive& DenseGapsRunner::writeSamplers(Archive &ar)
{
    ar << mASampler << mPSampler;
    return ar;
}

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
