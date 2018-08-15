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

GapsResult GapsRunner::run(bool printThreads)
{
    mStartTime = bpt_now();

    // calculate appropiate number of threads if compiled with openmp
    #ifdef __GAPS_OPENMP__
    if (mPrintMessages && printThreads)
    {
        unsigned availableThreads = omp_get_max_threads();
        mMaxThreads = gaps::min(availableThreads, mMaxThreads);
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
    GapsResult result(mStatistics);
    result.meanChiSq = mStatistics.meanChiSq(mASampler);
    return result;    
}

void GapsRunner::runOnePhase()
{
    for (; mCurrentIteration < mMaxIterations; ++mCurrentIteration)
    {
        createCheckpoint();

        #ifdef __GAPS_R_BUILD__
        Rcpp::checkUserInterrupt();
        #endif
        
        if (mPhase == 'C')
        {        
            float temp = static_cast<float>(2 * mCurrentIteration)
                / static_cast<float>(mMaxIterations);
            mASampler.setAnnealingTemp(gaps::min(1.f, temp));
            mPSampler.setAnnealingTemp(gaps::min(1.f, temp));
        }
    
        // number of updates per iteration is poisson 
        unsigned nA = mRng.poisson(gaps::max(mASampler.nAtoms(), 10));
        unsigned nP = mRng.poisson(gaps::max(mPSampler.nAtoms(), 10));
        updateSampler(nA, nP);

        if (mPhase == 'S')
        {
            mStatistics.update(mASampler, mPSampler);
        }

        displayStatus();
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
        GAPS_ASSERT(mASampler.internallyConsistent());
    }

    if (mFixedMatrix != 'P')
    {
        mNumUpdatesP += nP;
        mPSampler.update(nP, mMaxThreads);
        if (mFixedMatrix != 'A')
        {
            mASampler.sync(mPSampler);
        }
        GAPS_ASSERT(mPSampler.internallyConsistent());
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

double GapsRunner::estimatedPercentComplete() const
{
    double nIter = static_cast<double>(mCurrentIteration);
    double nAtomsA = static_cast<double>(mASampler.nAtoms());
    double nAtomsP = static_cast<double>(mPSampler.nAtoms());
    
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

void GapsRunner::displayStatus()
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
            mCurrentIteration + 1, mMaxIterations, mASampler.nAtoms(),
            mPSampler.nAtoms(), mASampler.chi2(), elapsedHours, elapsedMinutes,
            elapsedSeconds, totalHours, totalMinutes, totalSeconds);
        gaps_flush();
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
        GapsRng::save(ar);
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
