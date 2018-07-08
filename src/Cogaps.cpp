#include "GapsAssert.h"
#include "GibbsSampler.h"
#include "Matrix.h"
#include "Archive.h"
#include "InternalState.h"
#include "SIMD.h"

#include <Rcpp.h>
#include <ctime>
#include <fstream>
#include <cstdio>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// used to convert defined macro values into strings
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

// boost time helpers
namespace bpt = boost::posix_time;
#define bpt_now() bpt::microsec_clock::local_time()

// keeps track of when checkpoints are made
static bpt::ptime lastCheckpoint; 
static std::string checkpointFile;

// save the current internal state to a file
static void createCheckpoint(GapsInternalState &state)
{
    // create backup file
    std::rename(checkpointFile.c_str(), (checkpointFile + ".backup").c_str());

    // record starting time
    bpt::ptime start = bpt_now();

    // save state to file, write magic number at beginning
    std::string fname(checkpointFile);
    Archive ar(fname, ARCHIVE_WRITE);
    gaps::random::save(ar);
    ar << state.sampler.nFactor() << state.nEquil << state.nSample << state;
    ar.close();

    // display time it took to create checkpoint
    bpt::time_duration diff = bpt_now() - start;
    double elapsed = diff.total_milliseconds() / 1000.;
    Rprintf("created checkpoint in %.3f seconds\n", elapsed);

    // delete backup file
    std::remove((checkpointFile + ".backup").c_str());
}

static void updateSampler(GapsInternalState &state)
{
    state.nUpdatesA += state.nIterA;
    for (unsigned j = 0; j < state.nIterA; ++j)
    {
        state.sampler.update('A');
    }
    GAPS_ASSERT(state.sampler.internallyConsistent('A'));

    state.nUpdatesP += state.nIterP;
    for (unsigned j = 0; j < state.nIterP; ++j)
    {
        state.sampler.update('P');
    }
    GAPS_ASSERT(state.sampler.internallyConsistent('P'));
}

static void makeCheckpointIfNeeded(GapsInternalState &state)
{
    bpt::time_duration diff = bpt_now() - lastCheckpoint;
    int diff_sec = diff.total_milliseconds() / 1000;
    if (diff_sec > state.checkpointInterval && state.checkpointInterval > 0)
    {
        createCheckpoint(state);
        lastCheckpoint = bpt_now();
    }
}

static void storeSamplerInfo(GapsInternalState &state, Vector &atomsA,
Vector &atomsP, Vector &chi2)
{
    chi2[state.iter] = state.sampler.chi2();
    atomsA[state.iter] = state.sampler.totalNumAtoms('A');
    atomsP[state.iter] = state.sampler.totalNumAtoms('P');
    state.nIterA = gaps::random::poisson(std::max(atomsA[state.iter], 10.f));
    state.nIterP = gaps::random::poisson(std::max(atomsP[state.iter], 10.f));
}

static void displayStatus(GapsInternalState &state, const std::string &type,
unsigned nIterTotal)
{
    if ((state.iter + 1) % state.nOutputs == 0 && state.messages)
    {
        Rprintf("%s %d of %d, Atoms:%d(%d) Chi2 = %.2f\n", type.c_str(),
            state.iter + 1, nIterTotal, state.sampler.totalNumAtoms('A'),
            state.sampler.totalNumAtoms('P'), state.sampler.chi2());
    }
}

static void takeSnapshots(GapsInternalState &state)
{
    if (state.nSnapshots && !((state.iter+1)%(state.nSample/state.nSnapshots)))
    {
        state.snapshotsA.push_back(state.sampler.normedAMatrix().rMatrix());
        state.snapshotsP.push_back(state.sampler.normedPMatrix().rMatrix());
    }    
}

static void runBurnPhase(GapsInternalState &state)
{
    for (; state.iter < state.nEquil; ++state.iter)
    {
        makeCheckpointIfNeeded(state);
        float temp = ((float)state.iter + 2.f) / ((float)state.nEquil / 2.f);
        state.sampler.setAnnealingTemp(std::min(1.f,temp));
        updateSampler(state);
        displayStatus(state, "Equil: ", state.nEquil);
        storeSamplerInfo(state, state.nAtomsAEquil, state.nAtomsPEquil,
            state.chi2VecEquil);
    }
}

static void runCoolPhase(GapsInternalState &state)
{
    for (; state.iter < state.nEquilCool; ++state.iter)
    {
        makeCheckpointIfNeeded(state);
        updateSampler(state);
    }
}

static void runSampPhase(GapsInternalState &state)
{
    for (; state.iter < state.nSample; ++state.iter)
    {
        makeCheckpointIfNeeded(state);
        updateSampler(state);
        state.sampler.updateStatistics();
        if (state.nPumpSamples && !((state.iter + 1) % (state.nSample / state.nPumpSamples)))
        {
            state.sampler.updatePumpStatistics();
        }
        takeSnapshots(state);
        displayStatus(state, "Samp: ", state.nSample);
        storeSamplerInfo(state, state.nAtomsASample, state.nAtomsPSample,
            state.chi2VecSample);
    }
}

// execute the steps of the algorithm, return list to R
static Rcpp::List runCogaps(GapsInternalState &state)
{
    // reset the checkpoint timer
    lastCheckpoint = bpt_now();

    // cascade down the various phases of the algorithm
    // this allows for starting in the middle of the algorithm
    switch (state.phase)
    {
        case GAPS_BURN:
            runBurnPhase(state);
            state.iter = 0;
            state.phase = GAPS_COOL;

        case GAPS_COOL:
            runCoolPhase(state);
            state.iter = 0;
            state.phase = GAPS_SAMP;

        case GAPS_SAMP:
            runSampPhase(state);
    }

    // combine chi2 vectors
    Vector chi2Vec(state.chi2VecEquil);
    chi2Vec.concat(state.chi2VecSample);

    // print final chi-sq value
    float meanChiSq = state.sampler.meanChiSq();
    if (state.messages)
    {
        Rprintf("Chi-Squared of Mean: %.2f\n", meanChiSq);
    }

    return Rcpp::List::create(
        Rcpp::Named("Amean") = state.sampler.AMeanRMatrix(),
        Rcpp::Named("Asd") = state.sampler.AStdRMatrix(),
        Rcpp::Named("Pmean") = state.sampler.PMeanRMatrix(),
        Rcpp::Named("Psd") = state.sampler.PStdRMatrix(),
        Rcpp::Named("ASnapshots") = Rcpp::wrap(state.snapshotsA),
        Rcpp::Named("PSnapshots") = Rcpp::wrap(state.snapshotsP),
        Rcpp::Named("atomsAEquil") = state.nAtomsAEquil.rVec(),
        Rcpp::Named("atomsASamp") = state.nAtomsASample.rVec(),
        Rcpp::Named("atomsPEquil") = state.nAtomsPEquil.rVec(),
        Rcpp::Named("atomsPSamp") = state.nAtomsPSample.rVec(),
        Rcpp::Named("chiSqValues") = chi2Vec.rVec(),
        Rcpp::Named("randSeed") = state.seed,
        Rcpp::Named("numUpdates") = state.nUpdatesA + state.nUpdatesP,
        Rcpp::Named("meanChi2") = meanChiSq,
        Rcpp::Named("pumpStats") = state.sampler.pumpMatrix(),
        Rcpp::Named("meanPatternAssignment") = state.sampler.meanPattern()
    );
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, unsigned nEquil,
unsigned nEquilCool, unsigned nSample, unsigned nOutputs, unsigned nSnapshots,
float alphaA, float alphaP, float maxGibbmassA, float maxGibbmassP, int seed,
bool messages, bool singleCellRNASeq, char whichMatrixFixed,
const Rcpp::NumericMatrix &FP, unsigned checkpointInterval,
const std::string &cptFile, unsigned pumpThreshold, unsigned nPumpSamples)
{
    // set seed
    uint32_t seedUsed = static_cast<uint32_t>(seed);
    if (seed < 0)
    {
        bpt::ptime epoch(boost::gregorian::date(1970,1,1));
        bpt::time_duration diff = bpt_now() - epoch;
        seedUsed = static_cast<uint32_t>(diff.total_milliseconds() % 1000);
    }
    gaps::random::setSeed(seedUsed);

    // create internal state from parameters and run from there
    GapsInternalState state(D, S, nFactor, nEquil, nEquilCool, nSample,
        nOutputs, nSnapshots, alphaA, alphaP, maxGibbmassA, maxGibbmassP, seed,
        messages, singleCellRNASeq, whichMatrixFixed, FP, checkpointInterval,
        static_cast<PumpThreshold>(pumpThreshold), nPumpSamples);
    checkpointFile = cptFile;
    return runCogaps(state);
}

// TODO add checksum to verify D,S matrices
// [[Rcpp::export]]
Rcpp::List cogapsFromCheckpoint_cpp(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, const std::string &fileName,
const std::string &cptFile)
{   
    Archive ar(fileName, ARCHIVE_READ);
    gaps::random::load(ar);

    // read parameters needed to calculate the size of the internal state
    unsigned nFactor = 0, nEquil = 0, nSample = 0;
    ar >> nFactor >> nEquil >> nSample;
    
    // construct empty state of the correct size, populate from file
    GapsInternalState state(D, S, nFactor, nEquil, nSample);
    ar >> state;

    // run cogaps from this internal state
    checkpointFile = cptFile;
    return runCogaps(state);
}

// [[Rcpp::export]]
void displayBuildReport_cpp()
{
#if defined( __clang__ )
    Rcpp::Rcout << "Compiled with Clang\n";
#elif defined( __INTEL_COMPILER )
    Rcpp::Rcout << "Compiled with Intel ICC/ICPC\n";
#elif defined( __GNUC__ )
    Rcpp::Rcout << "Compiled with GCC v" << STR( __GNUC__ ) << "."
        << STR( __GNUC_MINOR__ ) << '\n';
#elif defined( _MSC_VER )
    Rcpp::Rcout << "Compiled with Microsoft Visual Studio\n";
#endif

#if defined( __GAPS_AVX__ )
    Rcpp::Rcout << "AVX enabled\n";
#elif defined( __GAPS_SSE__ )
    Rcpp::Rcout << "SSE enabled\n";
#else
    Rcpp::Rcout << "SIMD not enabled\n";
#endif
}
