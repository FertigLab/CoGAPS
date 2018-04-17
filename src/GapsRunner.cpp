#include "GapsRunner.h"

GapsRunner::GapsRunner(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
unsigned nFactor, unsigned nEquil, unsigned nCool, unsigned nSample,
unsigned nOutputs, unsigned nSnapshots, float alphaA, float alphaP,
float maxGibbsMassA, float maxGibbsMassP, uint32_t seed, bool messages,
bool singleCellRNASeq, unsigned cptInterval, const std::string &cptFile,
char whichMatrixFixed, const Rcpp::NumericMatrix &FP)
    :
mChiSqEquil(nEquil), mNumAAtomsEquil(nEquil), mNumPAtomsEquil(nEquil),
mChiSqSample(nSample), mNumAAtomsSample(nSample), mNumPAtomsSample(nSample),
mIterA(10), mIterP(10), mEquilIter(nEquil), mCoolIter(nCool),
mSampleIter(nSample), mNumPatterns(nFactor), mNumOutputs(nOutputs),
mPrintMessages(messages), mCurrentIter(0), mPhase(GAPS_BURN), mSeed(seed),
mCheckpointInterval(cptInterval), mCheckpointFile(cptFile),
mNumUpdatesA(0), mNumUpdatesP(0),
mASampler(D, S, nFactor, alphaA, maxGibbsMassA),
mPSampler(D, S, nFactor, alphaP, maxGibbsMassP),
mStatistics(D.nrow(), D.ncol(), nFactor)
{
    mASampler.sync(mPSampler);
    mPSampler.sync(mASampler);
    gaps::random::setSeed(seed);
}

GapsRunner::GapsRunner(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
unsigned nFactor, unsigned nEquil, unsigned nSample, const std::string &cptFile)
    :
mChiSqEquil(nEquil), mNumAAtomsEquil(nEquil), mNumPAtomsEquil(nEquil),
mChiSqSample(nSample), mNumAAtomsSample(nSample), mNumPAtomsSample(nSample),
mASampler(D, S, nFactor), mPSampler(D, S, nFactor),
mStatistics(D.nrow(), D.ncol(), nFactor)
{
    Archive ar(cptFile, ARCHIVE_READ);
    gaps::random::load(ar);

   //ar >> mChiSqEquil >> mNumAAtomsEquil >> mNumPAtomsEquil >> mChiSqSample
   //    >> mNumAAtomsSample >> mNumPAtomsSample >> mIterA >> mIterP
   //    >> mEquilIter >> mCoolIter >> mSampleIter >> mNumPatterns >> mNumOutputs
   //    >> mPrintMessages >> mCurrentIter >> mPhase >> mSeed
   //    >> mCheckpointInterval >> mCheckpointFile >> mNumUpdatesA
   //    >> mNumUpdatesP >> mASampler >> mPSampler >> mStatistics;

    mASampler.sync(mPSampler);
    mPSampler.sync(mASampler);
}

// execute the steps of the algorithm, return list to R
Rcpp::List GapsRunner::run()
{
    // reset the checkpoint timer
    mLastCheckpoint = bpt_now();

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
        Rcpp::Named("meanChi2") = meanChiSq
    );
}

void GapsRunner::runBurnPhase()
{
    for (; mCurrentIter < mEquilIter; ++mCurrentIter)
    {
        Rcpp::checkUserInterrupt();
        makeCheckpointIfNeeded();
        float temp = ((float)mCurrentIter + 2.f) / ((float)mEquilIter / 2.f);
        mASampler.setAnnealingTemp(std::min(1.f,temp));
        mPSampler.setAnnealingTemp(std::min(1.f,temp));
        updateSampler();
        displayStatus("Equil: ", mEquilIter);
        storeSamplerInfo(mNumAAtomsEquil, mNumPAtomsEquil, mChiSqEquil);
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
        displayStatus("Samp: ", mSampleIter);
        storeSamplerInfo(mNumAAtomsSample, mNumPAtomsSample, mChiSqSample);
    }
}

void GapsRunner::updateSampler()
{
    mNumUpdatesA += mIterA;
    mASampler.update(mIterA);
    mPSampler.sync(mASampler);

    mNumUpdatesP += mIterP;
    mPSampler.update(mIterP);
    mASampler.sync(mPSampler);
}

void GapsRunner::storeSamplerInfo(Vector &atomsA, Vector &atomsP, Vector &chi2)
{
    chi2[mCurrentIter] = mASampler.chi2();
    atomsA[mCurrentIter] = mASampler.nAtoms();
    atomsP[mCurrentIter] = mPSampler.nAtoms();
    mIterA = gaps::random::poisson(std::max(atomsA[mCurrentIter], 10.f));
    mIterP = gaps::random::poisson(std::max(atomsP[mCurrentIter], 10.f));
}

void GapsRunner::displayStatus(const std::string &type, unsigned nIterTotal)
{
    if ((mCurrentIter + 1) % mNumOutputs == 0 && mPrintMessages)
    {
        Rprintf("%s %d of %d, Atoms:%lu(%lu) Chi2 = %.2f\n", type.c_str(),
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
    //ar << mChiSqEquil << mNumAAtomsEquil << mNumPAtomsEquil << mChiSqSample
    //    << mNumAAtomsSample << mNumPAtomsSample << mIterA << mIterP
    //    << mEquilIter << mCoolIter << mSampleIter << mNumPatterns << mNumOutputs
    //    << mPrintMessages << mCurrentIter << mPhase << mSeed
    //    << mCheckpointInterval << mNumUpdatesA << mNumUpdatesP << mASampler
    //    << mPSampler << mStatistics;
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
    long sec = diff.total_milliseconds() / 1000;
    if (sec > mCheckpointInterval && mCheckpointInterval > 0)
    {
        createCheckpoint();
        mLastCheckpoint = bpt_now();
    }
}