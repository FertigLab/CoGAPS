#include "GapsDispatcher.h"

#include <string>

#ifdef __GAPS_OPENMP__
    #include <omp.h>
#endif

GapsDispatcher::GapsDispatcher() : mSeed(0), mNumPatterns(3),
    mMaxIterations(1000), mNumCoresPerSet(1), mPrintMessages(true),
    mCheckpointsCreated(0), mPhase('C'), mCheckpointInterval(0),
    mCheckpointOutFile("gaps_checkpoint.out"), mInitialized(false)
{}

GapsDispatcher::~GapsDispatcher()
{
    for (unsigned i = 0; i < mRunners.size(); ++i)
    {
        delete mRunners[i];
    }
}

void GapsDispatcher::setMaxIterations(unsigned n)
{
    mMaxIterations = n;
}

void GapsDispatcher::printMessages(bool print)
{
    GAPS_ASSERT(mInitialized);
    mPrintMessages = print;
    mRunners[0]->printMessages(print);
}

void GapsDispatcher::setOutputFrequency(unsigned n)
{
    GAPS_ASSERT(mInitialized);
    mRunners[0]->setOutputFrequency(n);
}

void GapsDispatcher::setSparsity(float alphaA, float alphaP, bool singleCell)
{
    GAPS_ASSERT(mInitialized);
    mRunners[0]->setSparsity(alphaA, alphaP, singleCell);
}

void GapsDispatcher::setMaxGibbsMass(float maxA, float maxP)
{
    GAPS_ASSERT(mInitialized);
    mRunners[0]->setMaxGibbsMass(maxA, maxP);
}

void GapsDispatcher::setAMatrix()
{
    GAPS_ASSERT(mInitialized);
}

void GapsDispatcher::setPMatrix()
{
    GAPS_ASSERT(mInitialized);
    mRunners[0]->
}

void GapsDispatcher::setFixedMatrix(char which)
{
    GAPS_ASSERT(mInitialized);
    mRunners[0]->setFixedMatrix(which, mat);
}

void GapsDispatcher::setNumCoresPerSet(unsigned n)
{
    mNumCoresPerSet = n;
}

void GapsDispatcher::setCheckpointInterval(unsigned n)
{
    mCheckpointInterval = n;
}

void GapsDispatcher::setCheckpointOutFile(const std::string &path)
{
    mCheckpointOutFile = path;
}

GapsResult GapsDispatcher::run()
{
    GAPS_ASSERT(mInitialized);
    GAPS_ASSERT(mPhase == 'C' || mPhase == 'S');

    // calculate appropiate number of cores if compiled with openmp
    #ifdef __GAPS_OPENMP__
    if (mPrintMessages)
    {
        unsigned availableCores = omp_get_max_threads();
        mNumCoresPerSet = std::min(availableCores, mNumCoresPerSet);
        gaps_printf("Running on %d out of %d available cores\n",
            mNumCoresPerSet, availableCores);
    }
    #endif

    // this switch allows for the algorithm to be interruptable
    mRunners[0]->startClock();
    switch (mPhase)
    {
        case 'C':
            if (mPrintMessages)
            {
                gaps_printf("-- Calibration Phase --\n");
            }
            runOneCycle(mMaxIterations);
            mPhase = 'S';

        case 'S':
            if (mPrintMessages)
            {
                gaps_printf("-- Sampling Phase --\n");
            }
            mRunners[0]->startSampling();
            runOneCycle(mMaxIterations);
            break;
    }

    // extract useful information from runners
    GapsResult result(mRunners[0]->nRow(), mRunners[0]->nCol(), mSeed);
    result.Amean = mRunners[0]->Amean();
    result.Pmean = mRunners[0]->Pmean();
    result.Asd = mRunners[0]->Asd();
    result.Psd = mRunners[0]->Psd();
    result.meanChiSq = mRunners[0]->meanChiSq();
    return result;
}

void GapsDispatcher::createCheckpoint() const
{
    Archive ar(mCheckpointOutFile, ARCHIVE_WRITE);

    gaps::random::save(ar);
    ar << mSeed << mNumPatterns << mMaxIterations << mPrintMessages <<
        mCheckpointsCreated << mPhase;

    ar << *mRunners[0];
}

void GapsDispatcher::runOneCycle(unsigned k)
{
    unsigned nCheckpoints = mCheckpointInterval > 0 ? k / mCheckpointInterval : 0;
    while (mCheckpointsCreated < nCheckpoints)
    {
        mRunners[0]->run(mCheckpointInterval, mNumCoresPerSet);
        createCheckpoint();
        ++mCheckpointsCreated;
    }
    mRunners[0]->run(k - mCheckpointInterval * mCheckpointsCreated, mNumCoresPerSet);
    mCheckpointsCreated = 0; // reset checkpoint count for next cycle
}
