#include "GapsAssert.h"
#include "GapsPrint.h"
#include "GapsRunner.h"
#include "math/Random.h"
#include "math/SIMD.h"

#ifdef __GAPS_R_BUILD__
#include <Rcpp.h>
#endif

#ifdef __GAPS_OPENMP__
#include <omp.h>
#endif

void GapsRunner::setFixedMatrix(char which, const Matrix &mat)
{
    mFixedMatrix = which;
    if (which == 'A')
    {
        mASampler.setMatrix(mat);
    }
    else if (which == 'P')
    {
        mPSampler.setMatrix(mat);
    }
}

void GapsRunner::recordSeed(uint32_t seed)
{
    mSeed = seed;
}

uint32_t GapsRunner::getSeed() const
{
    return mSeed;
}

void GapsRunner::setMaxIterations(unsigned nIterations)
{
    mMaxIterations = nIterations;
}

void GapsRunner::setSparsity(float alphaA, float alphaP, bool singleCell)
{
    mASampler.setSparsity(alphaA, singleCell);
    mPSampler.setSparsity(alphaP, singleCell);
}

void GapsRunner::setMaxGibbsMass(float maxA, float maxP)
{
    mASampler.setMaxGibbsMass(maxA);
    mPSampler.setMaxGibbsMass(maxP);
}

void GapsRunner::setMaxThreads(unsigned nThreads)
{
    mMaxThreads = nThreads;
}

void GapsRunner::setPrintMessages(bool print)
{
    mPrintMessages = print;
}

void GapsRunner::setOutputFrequency(unsigned n)
{
    mOutputFrequency = n;
}

void GapsRunner::setCheckpointOutFile(const std::string &file)
{
    mCheckpointOutFile = file;
}

void GapsRunner::setCheckpointInterval(unsigned interval)
{
    mCheckpointInterval = interval;
}

GapsResult GapsRunner::run()
{
#if 0
    mStartTime = bpt_now();

    // calculate appropiate number of threads if compiled with openmp
    #ifdef __GAPS_OPENMP__
    if (mPrintMessages)
    {
        unsigned availableThreads = omp_get_max_threads();
        mMaxThreads = std::min(availableThreads, mMaxThreads);
        gaps_printf("Running on %d out of %d available threads\n",
            mMaxThreads, availableThreads);
    }
    #endif

    // cascade through phases, allows algorithm to be resumed in either phase
    GAPS_ASSERT(mPhase == 'C' || mPhase == 'S');
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
#endif
    GapsResult result(mStatistics);
    result.meanChiSq = mStatistics.meanChiSq(mASampler);
    return result;    
}

void GapsRunner::runOnePhase()
{
    for (; mCurrentIteration < mMaxIterations; ++mCurrentIteration)
    {
        #ifdef __GAPS_R_BUILD__
        Rcpp::checkUserInterrupt();
        #endif
        
        if (mPhase == 'C')
        {        
            float temp = static_cast<float>(2 * mCurrentIteration)
                / static_cast<float>(mMaxIterations);
            mASampler.setAnnealingTemp(std::min(1.f, temp));
            mPSampler.setAnnealingTemp(std::min(1.f, temp));
        }
    
        // number of updates per iteration is poisson 
        unsigned nA = gaps::random::poisson(std::max(mASampler.nAtoms(), 10ul));
        unsigned nP = gaps::random::poisson(std::max(mPSampler.nAtoms(), 10ul));
        updateSampler(nA, nP);

        if (mPhase == 'S')
        {
            mStatistics.update(mASampler, mPSampler);
        }

        displayStatus();
        createCheckpoint();
    }
}

void GapsRunner::updateSampler(unsigned nA, unsigned nP)
{
    if (mFixedMatrix != 'A')
    {
        mNumUpdatesA += nA;
        mASampler.update(nA, mMaxThreads);
        if (mFixedMatrix != 'P')
        {
            mPSampler.sync(mASampler);
        }
    }

    if (mFixedMatrix != 'P')
    {
        mNumUpdatesP += nP;
        mPSampler.update(nP, mMaxThreads);
        if (mFixedMatrix != 'A')
        {
            mASampler.sync(mPSampler);
        }
    }

    GAPS_ASSERT(mASampler.internallyConsistent());
    GAPS_ASSERT(mPSampler.internallyConsistent());
}

void GapsRunner::displayStatus()
{
    if (mPrintMessages && mOutputFrequency > 0 && ((mCurrentIteration + 1) % mOutputFrequency) == 0)
    {
        bpt::time_duration diff = bpt_now() - mStartTime;
        unsigned elapsedSeconds = static_cast<unsigned>(diff.total_seconds());

        unsigned hours = elapsedSeconds / 3600;
        elapsedSeconds -= hours * 3600;
        unsigned minutes = elapsedSeconds / 60;
        elapsedSeconds -= minutes * 60;
        unsigned seconds = elapsedSeconds;

        gaps_printf("%d of %d, Atoms: %lu(%lu), ChiSq: %.0f, elapsed time: %02d:%02d:%02d\n",
            mCurrentIteration + 1, mMaxIterations, mASampler.nAtoms(),
            mPSampler.nAtoms(), mASampler.chi2(), hours, minutes, seconds);
    }
}

void GapsRunner::createCheckpoint()
{
    if (mCheckpointInterval > 0 && ((mCurrentIteration + 1) % mCheckpointInterval) == 0)
    {
        // create backup file
        std::rename(mCheckpointOutFile.c_str(), (mCheckpointOutFile + ".backup").c_str());
    
        // create checkpoint file
        Archive ar(mCheckpointOutFile, ARCHIVE_WRITE);
        
        gaps::random::save(ar);
        ar << mNumPatterns << mSeed << mASampler << mPSampler << mStatistics
            << mFixedMatrix << mMaxIterations << mPhase << mCurrentIteration
            << mNumUpdatesA << mNumUpdatesP;

        ar.close();

        // delete backup file
        std::remove((mCheckpointOutFile + ".backup").c_str());
    }
}

// assume random state has been loaded and nPatterns and seed have been read
Archive& operator>>(Archive &ar, GapsRunner &gr)
{
    ar >> gr.mNumPatterns >> gr.mSeed >> gr.mASampler >> gr.mPSampler
        >> gr.mStatistics >> gr.mFixedMatrix >> gr.mMaxIterations >> gr.mPhase
        >> gr.mCurrentIteration >> gr.mNumUpdatesA >> gr.mNumUpdatesP;
    return ar;
}