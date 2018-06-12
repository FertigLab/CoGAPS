#ifndef __COGAPS_GAPS_RUNNER_H__
#define __COGAPS_GAPS_RUNNER_H__

#include "Archive.h"
#include "GapsStatistics.h"
#include "GibbsSampler.h"

#include "data_structures/Matrix.h"
#include "data_structures/Vector.h"

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

// Manages all data and parameters for running CoGAPS; contains the high-level
// algorithm logic
class GapsRunner
{
private:
#ifdef GAPS_INTERNAL_TESTS
public:
#endif

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

    //unsigned mNumPatterns;
    unsigned mNumOutputs;
    bool mPrintMessages;

    unsigned mCurrentIter;
    GapsPhase mPhase;
    uint32_t mSeed;

    bpt::ptime mLastCheckpoint;
    int64_t mCheckpointInterval;
    std::string mCheckpointFile;

    unsigned mNumUpdatesA;
    unsigned mNumUpdatesP;
    
    AmplitudeGibbsSampler mASampler;
    PatternGibbsSampler mPSampler;
    GapsStatistics mStatistics;

    unsigned mNumCores;

    bpt::ptime mStartTime;

    void createCheckpoint();
    void makeCheckpointIfNeeded();
    void displayStatus(const std::string &type, unsigned nIterTotal);
    void storeSamplerInfo(Vector &atomsA, Vector &atomsP, Vector &chi2);
    void updateSampler();
    void runBurnPhase();
    void runCoolPhase();
    void runSampPhase();
    double estPercentComplete();

public:

    // construct from file name
    GapsRunner();

    // construct from R object
    GapsRunner(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
        unsigned nFactor, float alphaA, float alphaP, float maxGibbsMassA,
        float maxGibbsMassP, bool singleCellRNASeq
    
    // construct from checkpoint file
    GapsRunner(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
        unsigned nFactor, unsigned nEquil, unsigned nSample,
        const std::string &cptFile);

    void run();
    void halt();

    void startSampling();

    ColMatrix getAMatrix();
    RowMatrix getPMatrix();

    ColMatrix

    void setPMatrix(const RowMatrix &Pmaster, float weight);

    ColMatrix AMean() const;
    ColMatrix AStd() const;
    RowMatrix PMean() const;
    RowMatrix PStd() const;
};

#endif // __COGAPS_GAPS_RUNNER_H__