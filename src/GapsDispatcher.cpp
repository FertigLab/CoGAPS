#include "GapsDispatcher.h"

#include <string>

#ifdef __GAPS_OPENMP__
    #include <omp.h>
#endif

static std::vector< std::vector<unsigned> > sampleIndices(unsigned n, unsigned nSets)
{
    unsigned setSize = n / nSets;
    std::vector< std::vector<unsigned> > sampleIndices;
    std::vector<unsigned> toBeSampled;
    for (unsigned i = 0; i < n; ++i)
    {
        toBeSampled.push_back(i);
    }

    for (unsigned i = 0; i < (nSets - 1); ++i)
    {
        sampleIndices.push_back(gaps::random::sample(toBeSampled, setSize));
    }

    GAPS_ASSERT(!toBeSampled.empty());

    sampleIndices.push_back(toBeSampled);
    return sampleIndices;
}

void GapsDispatcher::runOneCycle(unsigned k)
{
    GAPS_ASSERT(mDataIsLoaded);
    mRunners[0]->run(k, mOutputFrequency, mPrintMessages, mNumCoresPerSet);
}

GapsResult GapsDispatcher::run()
{
    #ifdef __GAPS_OPENMP__
    if (mPrintMessages)
    {
        unsigned availableCores = omp_get_max_threads();
        gaps_printf("Running on %d out of %d available cores\n",
            mNumCoresPerSet, availableCores);
    }
    #endif

    GapsResult result(mRunners[0]->nRow(), mRunners[0]->nCol());
    GAPS_ASSERT(mDataIsLoaded);
    runOneCycle(mMaxIterations);
    mRunners[0]->startSampling();
    runOneCycle(mMaxIterations);

    result.Amean = mRunners[0]->AMean();
    result.Pmean = mRunners[0]->PMean();
    result.Asd = mRunners[0]->AStd();
    result.Psd = mRunners[0]->PStd();
    result.seed = mSeed;
    return result;
}

void GapsDispatcher::loadData(const RowMatrix &D)
{
    mRunners.push_back(new GapsRunner(D, mNumPatterns, mAlphaA,
        mAlphaP, mMaxGibbsMassA, mMaxGibbsMassP, mSingleCell));
    mDataIsLoaded = true;
}

void GapsDispatcher::setUncertainty(const RowMatrix &S)
{
    GAPS_ASSERT(mDataIsLoaded);
    mRunners[0]->setUncertainty(S);
}

void GapsDispatcher::loadData(const std::string &pathToData)
{
    FileParser parser(pathToData);
    mRunners.push_back(new GapsRunner(parser, mNumPatterns, mAlphaA,
        mAlphaP, mMaxGibbsMassA, mMaxGibbsMassP, mSingleCell));
    mDataIsLoaded = true;
}

void GapsDispatcher::setUncertainty(const std::string &pathToMatrix)
{
    GAPS_ASSERT(mDataIsLoaded);
    FileParser parser(pathToMatrix);
    mRunners[0]->setUncertainty(parser);
}