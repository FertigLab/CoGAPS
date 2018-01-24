#include "GibbsSampler.h"
#include "Matrix.h"
#include "Archive.h"
#include "InternalState.h"
#include "SIMD.h"

#include <Rcpp.h>
#include <ctime>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// no C++11 std::to_string
#include <sstream>
#define SSTR(x) static_cast<std::ostringstream&>( \
    (std::ostringstream() << std::dec << x)).str()

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#define ARCHIVE_MAGIC_NUM 0xCE45D32A

namespace bpt = boost::posix_time;
static bpt::ptime lastCheckpoint; // keep track of when checkpoints are made

// save the current internal state to a file
static void createCheckpoint(GapsInternalState &state)
{
    Rcpp::Rcout << "creating gaps checkpoint...";
    state.numCheckpoints++;

    // record starting time
    bpt::ptime start = bpt::microsec_clock::local_time();

    // save state to file, write magic number at beginning
    Archive ar("gaps_checkpoint_" + SSTR(state.numCheckpoints) + ".out",
        ARCHIVE_WRITE);
    ar << ARCHIVE_MAGIC_NUM;
    gaps::random::save(ar);
    ar << state.sampler.nFactor() << state.nEquil << state.nSample << state;
    ar.close();

    // display time it took to create checkpoint
    bpt::time_duration diff = bpt::microsec_clock::local_time() - start;
    double elapsed = diff.total_milliseconds() / 1000.;
    Rcpp::Rcout << " finished in " << elapsed << " seconds\n";
}

static void runGibbsSampler(GapsInternalState &state, unsigned nIterTotal,
Vector &chi2Vec, Vector &aAtomVec, Vector &pAtomVec)
{
    for (; state.iter < nIterTotal; ++state.iter)
    {
        bpt::ptime now = bpt::microsec_clock::local_time();
        bpt::time_duration diff = now - lastCheckpoint;
        if (diff.total_milliseconds() > state.checkpointInterval * 1000
        && state.checkpointInterval > 0)
        {
            createCheckpoint(state);
            lastCheckpoint = bpt::microsec_clock::local_time();
        }

        if (state.phase == GAPS_BURN)
        {
            state.sampler.setAnnealingTemp(std::min(1.f,
                ((float)state.iter + 2.f) / ((float)state.nEquil / 2.f)));
        }
        
        state.nUpdatesA += state.nIterA;
        state.nUpdatesP += state.nIterP;
        for (unsigned j = 0; j < state.nIterA; ++j)
        {
            state.sampler.update('A');
        }
        for (unsigned j = 0; j < state.nIterP; ++j)
        {
            state.sampler.update('P');
        }

        if (state.phase == GAPS_SAMP)
        {
            state.sampler.updateStatistics();
            if (state.nSnapshots && !((state.iter + 1) %
            (nIterTotal / state.nSnapshots)))
            {
                state.snapshotsA.push_back(state.sampler.getNormedMatrix('A'));
                state.snapshotsP.push_back(state.sampler.getNormedMatrix('P'));
            }
        }

        if (state.phase != GAPS_COOL)
        {
            float nAtomsA = state.sampler.totalNumAtoms('A');
            float nAtomsP = state.sampler.totalNumAtoms('P');
            aAtomVec[state.iter] = nAtomsA;
            pAtomVec[state.iter] = nAtomsP;
            state.nIterA = gaps::random::poisson(std::max(nAtomsA, 10.f));
            state.nIterP = gaps::random::poisson(std::max(nAtomsP, 10.f));

            if ((state.iter + 1) % state.nOutputs == 0 && state.messages)
            {
                std::string ph(state.phase == GAPS_BURN ? "Equil: " : "Samp: ");
                Rcpp::Rcout << ph << state.iter + 1 << " of " << nIterTotal
                    << ", Atoms:" << aAtomVec[state.iter] << "("
                    << pAtomVec[state.iter] << ") Chi2 = "
                    << state.sampler.chi2() << '\n';
            }
        }
    }
}

// execute the steps of the algorithm, return list to R
static Rcpp::List runCogaps(GapsInternalState &state)
{
    // reset the checkpoint timer
    lastCheckpoint = bpt::microsec_clock::local_time();

    // cascade down the various phases of the algorithm
    // this allows for starting in the middle of the algorithm
    Vector trash(1);
    switch (state.phase)
    {
        case GAPS_BURN:
            runGibbsSampler(state, state.nEquil, state.chi2VecEquil,
                state.nAtomsAEquil, state.nAtomsPEquil);
            state.iter = 0;
            state.phase = GAPS_COOL;
        case GAPS_COOL:
            runGibbsSampler(state, state.nEquilCool, trash, trash, trash);
            state.iter = 0;
            state.phase = GAPS_SAMP;
        case GAPS_SAMP:
            runGibbsSampler(state, state.nSample, state.chi2VecSample,
                state.nAtomsASample, state.nAtomsPSample);
    }

    // combine chi2 vectors
    Vector chi2Vec(state.chi2VecEquil);
    chi2Vec.concat(state.chi2VecSample);

    // print final chi-sq value
    float meanChiSq = state.sampler.meanChiSq();
    if (state.messages)
    {
        Rcpp::Rcout << "Chi-Squared of Mean: " << meanChiSq << '\n';
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
        Rcpp::Named("meanChi2") = meanChiSq
    );
}

// [[Rcpp::export]]
Rcpp::List cogaps_cpp(const Rcpp::NumericMatric &DMatrix,
const Rcpp::NumericMatrix &SMatrix, const Rcpp::NumericMatrix &fixedPatterns,
const Rcpp::S4 &params)
{
    // set seed
    uint32_t seedUsed = static_cast<uint32_t>(seed);
    if (seed < 0)
    {
        bpt::ptime epoch(boost::gregorian::date(1970,1,1)); 
        bpt::time_duration diff = bpt::microsec_clock::local_time() - epoch;
        seedUsed = static_cast<uint32_t>(diff.total_milliseconds() % 1000);
    }
    gaps::random::setSeed(seedUsed);

    // create internal state from parameters and run from there
    GapsInternalState state(DMatrix, SMatrix, fixedPatterns, params);
    return runCogaps(state);    
}

// [[Rcpp::export]]
Rcpp::List cogapsFromCheckpoint_cpp(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, const std::string &fileName)
{   
    // open file
    Archive ar(fileName, ARCHIVE_READ);

    // verify checkpoint file is valid by checking first 4 bytes
    // TODO add checksum to verify D,S matrices
    uint32_t magicNum = 0;
    ar >> magicNum;
    if (magicNum != ARCHIVE_MAGIC_NUM)
    {
        Rcpp::Rcout << "invalid checkpoint file" << std::endl;
        return Rcpp::List::create();
    }
    
    // seed random number generator
    gaps::random::load(ar);

    // read parameters needed to calculate the size of the internal state
    unsigned nFactor = 0, nEquil = 0, nSample = 0;
    ar >> nFactor >> nEquil >> nSample;
    
    // construct empty state of the correct size, populate from file
    GapsInternalState state(D, S, nFactor, nEquil, nSample);
    ar >> state;

    // run cogaps from this internal state
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
