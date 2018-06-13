#include "GapsAssert.h"
#include "GapsRunner.h"
#include "math/Random.h"
#include "math/SIMD.h"
#include <algorithm>

#ifdef __GAPS_OPENMP__
    #include <omp.h>
#endif

// create "nSets" vectors where each vector contains a vector of indices in the
// range [0,n)
static std::vector< std::vector<unsigned> > sampleIndices(unsigned n, unsigned nSets)
{
    unsigned setSize = n / nSets;
    std::vector< std::vector<unsigned> > sampleIndices;
    std::vector<unsigned> toBeSampled;
    for (unsigned i = 1; i < n; ++i)
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

GapsRunner::GapsRunner(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
unsigned nFactor, unsigned nEquil, unsigned nCool, unsigned nSample,
unsigned nOutputs, float alphaA, float alphaP,
float maxGibbsMassA, float maxGibbsMassP, uint32_t seed, bool messages,
bool singleCellRNASeq, unsigned cptInterval, const std::string &cptFile,
char whichMatrixFixed, const Rcpp::NumericMatrix &FP, unsigned nCores,
PumpThreshold pumpThreshold, unsigned nPumpSamples)
    :
mChiSqEquil(nEquil), mNumAAtomsEquil(nEquil), mNumPAtomsEquil(nEquil),
mChiSqSample(nSample), mNumAAtomsSample(nSample), mNumPAtomsSample(nSample),
mIterA(10), mIterP(10), mEquilIter(nEquil), mCoolIter(nCool),
mSampleIter(nSample), mNumOutputs(nOutputs),
mPrintMessages(messages), mCurrentIter(0), mPhase(GAPS_BURN), mSeed(seed),
mCheckpointInterval(cptInterval), mCheckpointFile(cptFile),
mNumUpdatesA(0), mNumUpdatesP(0),
mASampler(D, S, nFactor, alphaA, maxGibbsMassA, singleCellRNASeq),
mPSampler(D, S, nFactor, alphaP, maxGibbsMassP, singleCellRNASeq),
mStatistics(D.nrow(), D.ncol(), nFactor, pumpThreshold),
mNumCores(nCores), mFixedMatrix(whichMatrixFixed),
mNumPumpSamples(nPumpSamples)
{
    if (mFixedMatrix == 'A')
    {
        mASampler.setMatrix(ColMatrix(FP));
    }
    else if (mFixedMatrix == 'P')
    {
        mPSampler.setMatrix(RowMatrix(FP));
    }

    mASampler.sync(mPSampler);
    mPSampler.sync(mASampler);
    gaps::random::setSeed(seed);
}

GapsRunner::GapsRunner(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
unsigned nFactor, unsigned nEquil, unsigned nSample, const std::string &cptFile)
    :
mChiSqEquil(nEquil), mNumAAtomsEquil(nEquil), mNumPAtomsEquil(nEquil),
mChiSqSample(nSample), mNumAAtomsSample(nSample), mNumPAtomsSample(nSample),
mCheckpointFile(cptFile), mASampler(D, S, nFactor), mPSampler(D, S, nFactor),
mStatistics(D.nrow(), D.ncol(), nFactor)
{
    Archive ar(cptFile, ARCHIVE_READ);
    gaps::random::load(ar);
    unsigned storedPhase = 0;

    ar >> mChiSqEquil >> mNumAAtomsEquil >> mNumPAtomsEquil >> mChiSqSample
        >> mNumAAtomsSample >> mNumPAtomsSample >> mIterA >> mIterP
        >> mEquilIter >> mCoolIter >> mSampleIter >> mNumOutputs
        >> mPrintMessages >> mCurrentIter >> storedPhase >> mSeed
        >> mCheckpointInterval >> mNumUpdatesA
        >> mNumUpdatesP >> mASampler >> mPSampler >> mStatistics >> mNumCores;

    mPhase = static_cast<GapsPhase>(storedPhase);
    mASampler.sync(mPSampler);
    mPSampler.sync(mASampler);
}

// execute the steps of the algorithm, return list to R
Rcpp::List GapsRunner::run()
{
    #ifdef __GAPS_OPENMP__
        if (mPrintMessages)
        {
            unsigned availableCores = omp_get_max_threads();
            Rprintf("Running on %d out of %d available cores\n", mNumCores, availableCores);
        }
    #endif

    // reset the checkpoint timer
    mStartTime = bpt_now();
    mLastCheckpoint = mStartTime;

    // cascade down the various phases of the algorithm
    // this allows for starting in the middle of the algorithm
    switch (mPhase)
    {
        case GAPS_BURN:
            runBurnPhase();
            mCurrentIter = 0;
            mPhase = GAPS_COOL;
        case GAPS_COOL:
            runCoolPhase();
            mCurrentIter = 0;
            mPhase = GAPS_SAMP;
        case GAPS_SAMP:
            runSampPhase();
    }

    // combine chi2 vectors
    Vector chi2Vec(mChiSqEquil);
    chi2Vec.concat(mChiSqSample);

    // print final chi-sq value
    float meanChiSq = mStatistics.meanChiSq(mASampler);
    if (mPrintMessages)
    {
        Rprintf("Chi-Squared of Mean: %.2f\n", meanChiSq);
    }

    return Rcpp::List::create(
        Rcpp::Named("Amean") = mStatistics.AMean(),
        Rcpp::Named("Asd") = mStatistics.AStd(),
        Rcpp::Named("Pmean") = mStatistics.PMean(),
        Rcpp::Named("Psd") = mStatistics.PStd(),
        Rcpp::Named("atomsAEquil") = mNumAAtomsEquil.rVec(),
        Rcpp::Named("atomsASamp") = mNumAAtomsSample.rVec(),
        Rcpp::Named("atomsPEquil") = mNumPAtomsEquil.rVec(),
        Rcpp::Named("atomsPSamp") = mNumPAtomsSample.rVec(),
        Rcpp::Named("chiSqValues") = chi2Vec.rVec(),
        Rcpp::Named("randSeed") = mSeed,
        Rcpp::Named("numUpdates") = mNumUpdatesA + mNumUpdatesP,
        Rcpp::Named("meanChi2") = meanChiSq,
        Rcpp::Named("AAvgQueue") = mASampler.getAvgQueue(),
        Rcpp::Named("PAvgQueue") = mPSampler.getAvgQueue(),
        Rcpp::Named("pumpStats") = mStatistics.pumpMatrix(),
        Rcpp::Named("meanPatternAssignment") = mStatistics.meanPattern()
    );
}

void GapsRunner::runBurnPhase()
{
    for (; mCurrentIter < mEquilIter; ++mCurrentIter)
    {
        Rcpp::checkUserInterrupt();
        makeCheckpointIfNeeded();
        float temp = (static_cast<float>(mCurrentIter) + 2.f) /
            (static_cast<float>(mEquilIter) / 2.f);
        mASampler.setAnnealingTemp(std::min(1.f,temp));
        mPSampler.setAnnealingTemp(std::min(1.f,temp));
        updateSampler();
        storeSamplerInfo(mNumAAtomsEquil, mNumPAtomsEquil, mChiSqEquil);
        displayStatus("Equil: ", mEquilIter);
    }
}

void GapsRunner::runCoolPhase()
{
    for (; mCurrentIter < mCoolIter; ++mCurrentIter)
    {
        Rcpp::checkUserInterrupt();
        makeCheckpointIfNeeded();
        updateSampler();
    }
}

void GapsRunner::runSampPhase()
{
    for (; mCurrentIter < mSampleIter; ++mCurrentIter)
    {
        Rcpp::checkUserInterrupt();
        makeCheckpointIfNeeded();
        updateSampler();
        mStatistics.update(mASampler, mPSampler);
        if (mNumPumpSamples != 0 && ((mCurrentIter + 1) % (mSampleIter / mNumPumpSamples)) == 0)
        {
            mStatistics.updatePump(mASampler, mPSampler);
        }
        storeSamplerInfo(mNumAAtomsSample, mNumPAtomsSample, mChiSqSample);
        displayStatus("Samp: ", mSampleIter);
    }
}

// sum coef * log(i) for i = 1 to total, fit coef from number of atoms
// approximates sum of number of atoms
// this should be proportional to total number of updates
static double estNumUpdates(double current, double total, float nAtoms)
{
    double coef = nAtoms / std::log(current);
    return coef * std::log(std::sqrt(2 * total * gaps::algo::pi)) +
        total * coef * std::log(total) - total * coef;
}

double GapsRunner::estPercentComplete()
{
    unsigned nIter = 0;
    float nAtomsA = 0.f, nAtomsP = 0.f;
    if (mPhase == GAPS_BURN)
    {
        nIter = mCurrentIter;
        nAtomsA = mNumAAtomsEquil[mCurrentIter];
        nAtomsP = mNumPAtomsEquil[mCurrentIter];
    }
    else if (mPhase == GAPS_SAMP)
    {
        nIter = mEquilIter + mCoolIter + mCurrentIter;
        nAtomsA = mNumAAtomsSample[mCurrentIter];
        nAtomsP = mNumPAtomsSample[mCurrentIter];
    }

    unsigned totalIter = mEquilIter + mCoolIter + mSampleIter;

    return (estNumUpdates(nIter, nIter, nAtomsA) +
        estNumUpdates(nIter, nIter, nAtomsP)) /
        (estNumUpdates(nIter, totalIter, nAtomsA) +
        estNumUpdates(nIter, totalIter, nAtomsP));
}

void GapsRunner::updateSampler()
{
    if (mFixedMatrix != 'A')
    {
        mNumUpdatesA += mIterA;
        mASampler.update(mIterA, mNumCores);
        mPSampler.sync(mASampler);
    }

    if (mFixedMatrix != 'P')
    {
        mNumUpdatesP += mIterP;
        mPSampler.update(mIterP, mNumCores);
        mASampler.sync(mPSampler);
    }
}

void GapsRunner::storeSamplerInfo(Vector &atomsA, Vector &atomsP, Vector &chi2)
{
    chi2[mCurrentIter] = mASampler.chi2();
    atomsA[mCurrentIter] = mASampler.nAtoms();
    atomsP[mCurrentIter] = mPSampler.nAtoms();
    mIterA = gaps::random::poisson(std::max(atomsA[mCurrentIter], 10.f));
    mIterP = gaps::random::poisson(std::max(atomsP[mCurrentIter], 10.f));
}

static void printTime(const std::string &message, unsigned totalSeconds)
{
    int hours = totalSeconds / 3600;
    totalSeconds -= hours * 3600;
    int minutes = totalSeconds / 60;
    totalSeconds -= minutes * 60;
    int seconds = totalSeconds;
    Rprintf("%s: %02d:%02d:%02d\n", message.c_str(), hours, minutes, seconds);
}

// need to call storeSamplerInfo before calling this function
void GapsRunner::displayStatus(const std::string &type, unsigned nIterTotal)
{
    if (mNumOutputs > 0 && (mCurrentIter + 1) % mNumOutputs == 0 && mPrintMessages)
    {
        // compute time
        bpt::time_duration diff = bpt_now() - mStartTime;
        double elapsed = diff.total_seconds();
        double estTotal = elapsed / estPercentComplete();

        // print time messages
        printTime("Elapsed Time", elapsed);
        printTime("Estimated Total Time", estTotal);
        printTime("Estimated Remaining Time", estTotal - elapsed);

        // print algo status
        Rprintf("%s %d of %d, Atoms:%lu(%lu) Chi2 = %.2f\n\n", type.c_str(),
            mCurrentIter + 1, nIterTotal, mASampler.nAtoms(),
            mPSampler.nAtoms(), mASampler.chi2());
    }
}

// save the current internal state to a file
void GapsRunner::createCheckpoint()
{
    // create backup file
    std::rename(mCheckpointFile.c_str(), (mCheckpointFile + ".backup").c_str());

    // record starting time
    bpt::ptime start = bpt_now();

    // save state to file, write magic number at beginning
    Archive ar(mCheckpointFile, ARCHIVE_WRITE);
    gaps::random::save(ar);
    ar << mChiSqEquil << mNumAAtomsEquil << mNumPAtomsEquil << mChiSqSample
        << mNumAAtomsSample << mNumPAtomsSample << mIterA << mIterP
        << mEquilIter << mCoolIter << mSampleIter << mNumOutputs
        << mPrintMessages << mCurrentIter << static_cast<unsigned>(mPhase)
        << mSeed 
        << mCheckpointInterval << mNumUpdatesA
        << mNumUpdatesP << mASampler << mPSampler << mStatistics << mNumCores;
    ar.close();

    // display time it took to create checkpoint
    bpt::time_duration diff = bpt_now() - start;
    double elapsed = diff.total_milliseconds() / 1000.;
    Rprintf("created checkpoint in %.3f seconds\n", elapsed);

    // delete backup file
    std::remove((mCheckpointFile + ".backup").c_str());
}

void GapsRunner::makeCheckpointIfNeeded()
{
    bpt::time_duration diff = bpt_now() - mLastCheckpoint;
    int64_t sec = diff.total_milliseconds() / 1000;
    if (sec > mCheckpointInterval && mCheckpointInterval > 0)
    {
        createCheckpoint();
        mLastCheckpoint = bpt_now();
    }
}
