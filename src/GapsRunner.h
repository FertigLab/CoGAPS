#ifndef __GAPS_GAPS_RUNNER_H__
#define __GAPS_GAPS_RUNNER_H__

#include "Archive.h"
#include "Matrix.h"
#include "GibbsSampler.h"
#include "GapsStatistics.h"

#include <Rcpp.h>

// boost time helpers
#include <boost/date_time/posix_time/posix_time.hpp>
namespace bpt = boost::posix_time;
#define bpt_now() bpt::microsec_clock::local_time()

enum GapsPhase
{
    GAPS_BURN,
    GAPS_COOL,
    GAPS_SAMP
};

class GapsRunner
{
private:

    Vector mChiSqEquil;
    Vector mNumAAtomsEquil;
    Vector mNumPAtomsEquil;

    Vector mChiSqSample;
    Vector mNumAAtomsSample;
    Vector mNumPAtomsSample;

    unsigned mIterA;
    unsigned mIterP;
    
    unsigned mEquilIter;
    unsigned mCoolIter;
    unsigned mSampleIter;

    unsigned mNumPatterns;
    unsigned mNumOutputs;
    bool mPrintMessages;

    unsigned mCurrentIter;
    GapsPhase mPhase;
    uint32_t mSeed;

    bpt::ptime mLastCheckpoint;
    long mCheckpointInterval;
    std::string mCheckpointFile;

    unsigned mNumUpdatesA;
    unsigned mNumUpdatesP;
    
    AmplitudeGibbsSampler mASampler;
    PatternGibbsSampler mPSampler;
    GapsStatistics mStatistics;

    void createCheckpoint();
    void makeCheckpointIfNeeded();
    void displayStatus(const std::string &type, unsigned nIterTotal);
    void storeSamplerInfo(Vector &atomA, Vector &atomsP, Vector &chi2);
    void updateSampler();
    void runBurnPhase();
    void runCoolPhase();
    void runSampPhase();

public:

    GapsRunner(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
        unsigned nFactor, unsigned nEquil, unsigned nCool, unsigned nSample,
        unsigned nOutputs, unsigned nSnapshots, float alphaA, float alphaP,
        float maxGibbsMassA, float maxGibbsMassP, uint32_t seed, bool messages,
        bool singleCellRNASeq, unsigned cptInterval, const std::string &cptFile,
        char whichMatrixFixed, const Rcpp::NumericMatrix &FP);

    GapsRunner(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
        unsigned nFactor, unsigned nEquil, unsigned nSample,
        const std::string &cptFile);

    Rcpp::List run();
};

#endif