#include "GapsRunner.h"
#include "GapsResult.h"
#include "GapsParameters.h"
#include "GapsStatistics.h"
#include "math/Random.h"
#include "utils/Archive.h"
#include "utils/GlobalConfig.h"
#include "gibbs_sampler/AsynchronousGibbsSampler.h"
#include "gibbs_sampler/SingleThreadedGibbsSampler.h"
#include "gibbs_sampler/DenseNormalModel.h"
#include "gibbs_sampler/SparseNormalModel.h"

#ifdef __GAPS_R_BUILD__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Rcpp.h>
#pragma GCC diagnostic pop
#endif

#ifdef __GAPS_OPENMP__
#include <omp.h>
#endif

// boost time helpers
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/date_time/posix_time/posix_time.hpp>
#pragma GCC diagnostic pop
namespace bpt = boost::posix_time;
#define bpt_now() bpt::microsec_clock::local_time()

// for conditionally printing status messages
#define GAPS_MESSAGE(b, m) \
    do { \
        if (b) { \
            gaps_printf(m); \
        } \
    } while(0)

// for converting seconds to h:m:s
struct GapsTime
{
    unsigned hours;
    unsigned minutes;
    unsigned seconds;
    explicit GapsTime(unsigned totalSeconds)
    {
        hours = totalSeconds / 3600;
        totalSeconds -= hours * 3600;
        minutes = totalSeconds / 60;
        totalSeconds -= minutes * 60;
        seconds = totalSeconds;
    }
};

// forward declaration
template <class Sampler, class DataType>
static GapsResult runCoGAPSAlgorithm(const DataType &data, GapsParameters &params,
    const DataType &uncertainty, GapsRandomState *randState);

////////////////////////////////////////////////////////////////////////////////

template <class DataModel, class DataType>
static GapsResult chooseSampler(const DataType &data, GapsParameters &params,
const DataType &uncertainty, GapsRandomState *randState)
{
    if (params.asynchronousUpdates)
    {
        GAPS_MESSAGE(params.printMessages, "Sampler Type: Asynchronous\n");
        return runCoGAPSAlgorithm< AsynchronousGibbsSampler<DataModel> >(data,
            params, uncertainty, randState);
    }
    GAPS_MESSAGE(params.printMessages, "Sampler Type: Sequential\n");
    return runCoGAPSAlgorithm< SingleThreadedGibbsSampler<DataModel> >(data,
        params, uncertainty, randState);
}

template <class DataType>
static GapsResult chooseDataModel(const DataType &data, GapsParameters &params,
const DataType &uncertainty, GapsRandomState *randState)
{
    if (params.useSparseOptimization)
    {
        GAPS_MESSAGE(params.printMessages, "Data Model: Sparse, Normal\n");
        return chooseSampler<SparseNormalModel>(data, params, uncertainty, randState);
    }
    GAPS_MESSAGE(params.printMessages, "Data Model: Dense, Normal\n");        
    return chooseSampler<DenseNormalModel>(data, params, uncertainty, randState);
}

// helper function, this dispatches the correct run function depending
// on the type of GibbsSampler and DataModel needed for the given parameters
template <class DataType>
static GapsResult run_helper(const DataType &data, GapsParameters &params,
const DataType &uncertainty, GapsRandomState *randState)
{
    // fetch parameters from checkpoint - some are used in initialization
    if (params.useCheckPoint)
    {
        Archive ar(params.checkpointFile, ARCHIVE_READ);
        ar >> params;
        ar >> *randState;
    }
    return chooseDataModel(data, params, uncertainty, randState);
}

// these two functions are the top-level functions exposed to the C++
// code that is being wrapped by any given language

GapsResult gaps::run(const Matrix &data, GapsParameters &params,
const Matrix &uncertainty, GapsRandomState *randState)
{
    return run_helper(data, params, uncertainty, randState);
}

GapsResult gaps::run(const std::string &data, GapsParameters &params,
const std::string &uncertainty, GapsRandomState *randState)
{
    return run_helper(data, params, uncertainty, randState);
}

////////////////////////////////////////////////////////////////////////////////

// sum coef * log(i) for i = 1 to total, fit coef from number of atoms
// approximates sum of number of atoms (stirling approx to factorial)
// this should be proportional to total number of updates
static double estimatedNumUpdates(double current, double total, float nAtoms)
{
    double coef = nAtoms / std::log(current);
    return coef * std::log(std::sqrt(2 * total * gaps::pi)) +
        total * coef * std::log(total) - total * coef;
}

template <class Sampler>
static double estimatedPercentComplete(const GapsParameters &params,
const Sampler &ASampler, const Sampler &PSampler, GapsAlgorithmPhase phase, unsigned iter)
{
    double nIter = static_cast<double>(iter);
    double nAtomsA = static_cast<double>(ASampler.nAtoms());
    double nAtomsP = static_cast<double>(PSampler.nAtoms());
    
    if (phase == GAPS_SAMPLING_PHASE)
    {
        nIter += params.nIterations;
    }

    double totalIter = 2.0 * static_cast<double>(params.nIterations);
    double estimatedCompleted = estimatedNumUpdates(nIter, nIter, nAtomsA) + 
        estimatedNumUpdates(nIter, nIter, nAtomsP);
    double estimatedTotal = estimatedNumUpdates(nIter, totalIter, nAtomsA) + 
        estimatedNumUpdates(nIter, totalIter, nAtomsP);
    return estimatedCompleted / estimatedTotal;
}

template <class Sampler>
static void displayStatus(const GapsParameters &params,
const Sampler &ASampler, const Sampler &PSampler, bpt::ptime startTime,
GapsAlgorithmPhase phase, unsigned iter, GapsStatistics &stats)
{
    if (params.outputFrequency > 0 && ((iter + 1) % params.outputFrequency) == 0)
    {
        float cs = PSampler.chiSq();
        unsigned nA = ASampler.nAtoms();
        unsigned nP = PSampler.nAtoms();
        stats.addChiSq(cs);
        stats.addAtomCount(nA, nP);
        if (params.printMessages)
        {
            bpt::time_duration diff = bpt_now() - startTime;
            double perComplete = estimatedPercentComplete(params, ASampler,
                PSampler, phase, iter);

            GapsTime elapsedTime(static_cast<unsigned>(diff.total_seconds()));
            GapsTime totalTime(static_cast<unsigned>(diff.total_seconds() / perComplete));

            gaps_printf("%d of %d, Atoms: %d(A), %d(P), ChiSq: %.0f, Time: %02d:%02d:%02d / %02d:%02d:%02d\n",
                iter + 1, params.nIterations, nA, nP, cs, elapsedTime.hours,
                elapsedTime.minutes, elapsedTime.seconds, totalTime.hours,
                totalTime.minutes, totalTime.seconds);
            gaps_flush();
        }
    }
}

template <class Sampler>
static void updateSampler(const GapsParameters &params, Sampler &ASampler,
Sampler &PSampler, unsigned nA, unsigned nP)
{
    if (params.whichMatrixFixed != 'A')
    {
        ASampler.update(nA, params.maxThreads);
        if (params.whichMatrixFixed != 'P')
        {
            PSampler.sync(ASampler, params.maxThreads);
        }
    }

    if (params.whichMatrixFixed != 'P')
    {
        PSampler.update(nP, params.maxThreads);
        if (params.whichMatrixFixed != 'A')
        {
            ASampler.sync(PSampler, params.maxThreads);
        }
    }
}

template <class Sampler>
static void createCheckpoint(const GapsParameters &params,
Sampler &ASampler, Sampler &PSampler, const GapsRandomState *randState,
const GapsStatistics &stats, const GapsRng &rng, GapsAlgorithmPhase phase, unsigned iter)
{
    if (params.checkpointInterval > 0 && ((iter + 1) % params.checkpointInterval) == 0
    && !params.subsetData)
    {
        // create backup file
        std::rename(params.checkpointOutFile.c_str(),
            (params.checkpointOutFile + ".backup").c_str());
    
        // create checkpoint file
        Archive ar(params.checkpointOutFile, ARCHIVE_WRITE);
        ar << params;
        ar << *randState;
        ar << ASampler << PSampler << stats << static_cast<int>(phase) << iter << rng;
        
        // delete backup file
        std::remove((params.checkpointOutFile + ".backup").c_str());

        // running the extra initialization here allows for consistency with runs
        // started from a checkpoint. This initialization phase will be run first 
        // thing once a checkpoint is loaded since large matrices which aren't stored
        // need to be initialized. By running it here we make sure that the algorithm
        // is in the same state it will be when started from a checkpoint
        ASampler.extraInitialization();
        PSampler.extraInitialization();
    }
}

template <class Sampler>
static void processCheckpoint(GapsParameters &params, Sampler &ASampler,
Sampler &PSampler, GapsRandomState *randState, GapsStatistics &stats,
GapsRng &rng, GapsAlgorithmPhase &phase, unsigned &currentIter)
{    
    // check if running from checkpoint, get all saved data
    if (params.useCheckPoint)
    {
        int iPhase;
        Archive ar(params.checkpointFile, ARCHIVE_READ);
        ar >> params;
        ar >> *randState;
        ar >> ASampler >> PSampler >> stats >> iPhase >> currentIter >> rng;
        phase = static_cast<GapsAlgorithmPhase>(iPhase);
    }
}

template <class Sampler>
static uint64_t runOnePhase(const GapsParameters &params, Sampler &ASampler,
Sampler &PSampler, GapsStatistics &stats, const GapsRandomState *randState,
GapsRng &rng, bpt::ptime startTime, GapsAlgorithmPhase phase, unsigned &currentIter)
{
    uint64_t totalUpdates = 0;
    for (; currentIter < params.nIterations; ++currentIter)
    {
        gaps_check_interrupt();
        createCheckpoint(params, ASampler, PSampler, randState, stats,
            rng, phase, currentIter);
        
        // set annealing temperature in Equilibration phase
        if (phase == GAPS_EQUILIBRATION_PHASE)
        {        
            float temp = static_cast<float>(2 * currentIter)
                / static_cast<float>(params.nIterations);
            ASampler.setAnnealingTemp(gaps::min(1.f, temp));
            PSampler.setAnnealingTemp(gaps::min(1.f, temp));
        }
    
        // number of updates per iteration is poisson 
        unsigned nA = rng.poisson(gaps::max(ASampler.nAtoms(), 10));
        unsigned nP = rng.poisson(gaps::max(PSampler.nAtoms(), 10));
        updateSampler(params, ASampler, PSampler, nA, nP);
        totalUpdates += nA + nP;

        if (phase == GAPS_SAMPLING_PHASE)
        {
            stats.update(ASampler, PSampler);
            if (params.takePumpSamples)
            {
                stats.updatePump(ASampler);
            }
        }
        if (params.snapshotPhase == phase || params.snapshotPhase == GAPS_ALL_PHASES)
        {
            if (params.snapshotFrequency > 0 && ((currentIter + 1) % params.snapshotFrequency) == 0)
            {
                stats.takeSnapshot(phase, ASampler, PSampler);
            }
        }
        displayStatus(params, ASampler, PSampler, startTime, phase,
            currentIter, stats);
    }
    return totalUpdates;
}

template <class Sampler>
static void processFixedMatrix(const GapsParameters &params, Sampler &ASampler,
Sampler &PSampler)
{
    // check if we're fixing a matrix
    if (params.useFixedPatterns)
    {
        GAPS_ASSERT(params.fixedPatterns.nCol() == params.nPatterns);
        switch (params.whichMatrixFixed)
        {
            case 'A' :
                GAPS_ASSERT(params.fixedPatterns.nRow() == params.nGenes);
                ASampler.setMatrix(params.fixedPatterns);
                break;
            case 'P' :
                GAPS_ASSERT(params.fixedPatterns.nRow() == params.nSamples);
                PSampler.setMatrix(params.fixedPatterns);
                break;
            default: break; // 'N' for none
        }
    }
}

static void calculateNumberOfThreads(GapsParameters params)
{
    // calculate appropiate number of threads if compiled with openmp
    #ifdef __GAPS_OPENMP__
    if (params.printMessages && params.printThreadUsage)
    {
        unsigned availableThreads = omp_get_max_threads();
        params.maxThreads = gaps::min(availableThreads, params.maxThreads);
        gaps_printf("Running on %d out of %d available threads\n",
            params.maxThreads, availableThreads);
    }
    #endif
}

template <class Sampler, class DataType>
static void processUncertainty(const GapsParameters params, Sampler &ASampler,
Sampler &PSampler, const DataType &uncertainty)
{
    // read in the uncertainty matrix if one is provided
    if (!uncertainty.empty())
    {
        ASampler.setUncertainty(uncertainty, !params.transposeData,
            !params.subsetGenes, params);
        PSampler.setUncertainty(uncertainty, params.transposeData,
            params.subsetGenes, params);
    }
}

// here is the CoGAPS algorithm
template <class Sampler, class DataType>
static GapsResult runCoGAPSAlgorithm(const DataType &data, GapsParameters &params,
const DataType &uncertainty, GapsRandomState *randState)
{
    // check if running in debug mode
    #ifdef GAPS_DEBUG
    gaps_printf("Running in debug mode\n");
    #endif

    // load data into gibbs samplers
    // we transpose the data in the A sampler so that the update step
    // is symmetrical for each sampler, this simplifies the code 
    // within the sampler, note the subsetting genes/samples flag must be
    // flipped if we are flipping the transpose flag
    // note: there are excessive amounts of check interrupt statements here,
    // this is so that if the user accidentally runs cogaps on an extremely
    // large file they can exit during the loading phase rather than killing
    // the process or waiting for it to finish
    GAPS_MESSAGE(params.printMessages, "Loading Data...");
    bpt::ptime readStart = bpt_now();
    gaps_check_interrupt();
    Sampler ASampler(data, !params.transposeData, !params.subsetGenes,
        params.alphaA, params.maxGibbsMassA, params, randState);
    gaps_check_interrupt();
    Sampler PSampler(data, params.transposeData, params.subsetGenes,
        params.alphaP, params.maxGibbsMassP, params, randState);
    gaps_check_interrupt();
    processUncertainty(params, ASampler, PSampler, uncertainty);
    gaps_check_interrupt();
    processFixedMatrix(params, ASampler, PSampler);
    gaps_check_interrupt();

    // elapsed time for reading data
    bpt::time_duration readDiff = bpt_now() - readStart;
    GapsTime elapsed(static_cast<unsigned>(readDiff.total_seconds()));
    if (params.printMessages)
    {
        gaps_printf("Done! (%02d:%02d:%02d)\n", elapsed.hours, elapsed.minutes,
            elapsed.seconds);
    }

    // check if data is sparse and sparseOptimization is not enabled
    if (params.printMessages && !params.useSparseOptimization && ASampler.dataSparsity() > 0.80f)
    {
        gaps_printf("\nWarning: data is more than 80%% sparse and sparseOptimization is not enabled\n");
    }

    // if we are running distributed, each worker needs to print when it's started
    if (params.runningDistributed)
    {
        gaps_printf("    worker %d is starting!\n", params.workerID);
        gaps_flush();
    }

    // these variables will get overwritten by checkpoint if provided
    GapsStatistics stats(params.nGenes, params.nSamples, params.nPatterns);
    GapsRng rng(randState);
    GapsAlgorithmPhase phase(GAPS_EQUILIBRATION_PHASE);
    unsigned currentIter = 0;
    processCheckpoint(params, ASampler, PSampler, randState, stats, rng, phase, currentIter);
    calculateNumberOfThreads(params);

    // sync samplers and run any additional initialization needed
    ASampler.sync(PSampler);
    PSampler.sync(ASampler);
    ASampler.extraInitialization();
    PSampler.extraInitialization();

    // record start time
    bpt::ptime startTime = bpt_now();

    // fallthrough through phases, allows algorithm to be resumed in either phase
    GAPS_ASSERT(phase == GAPS_EQUILIBRATION_PHASE || phase == GAPS_SAMPLING_PHASE);
    uint64_t totalUpdates = 0;
    switch (phase)
    {
        case GAPS_EQUILIBRATION_PHASE:
            GAPS_MESSAGE(params.printMessages, "-- Equilibration Phase --\n");
            totalUpdates += runOnePhase(params, ASampler, PSampler, stats, randState,
                rng, startTime, phase, currentIter);
            phase = GAPS_SAMPLING_PHASE;
            currentIter = 0;
        // fall through
        case GAPS_SAMPLING_PHASE:
            GAPS_MESSAGE(params.printMessages, "-- Sampling Phase --\n");
            totalUpdates += runOnePhase(params, ASampler, PSampler, stats, randState,
                rng, startTime, phase, currentIter);
    }
    
    // get result
    GapsResult result(stats);
    result.totalRunningTime = static_cast<unsigned>((bpt_now() - startTime).total_seconds());
    result.meanChiSq = stats.meanChiSq(PSampler);
    result.averageQueueLengthA = ASampler.getAverageQueueLength();
    result.averageQueueLengthP = PSampler.getAverageQueueLength();
    result.totalUpdates = totalUpdates;

    // handle pump statistics
    if (params.takePumpSamples)
    {
        result.pumpMatrix = stats.pumpMatrix();
        result.meanPatternAssignment = stats.meanPattern();
    }

    // if we are running distributed, each worker needs to print when it's done
    if (params.runningDistributed)
    {
        GapsTime elapsed(static_cast<unsigned>((bpt_now() - startTime).total_seconds()));
        gaps_printf("    worker %d is finished! Time: %02d:%02d:%02d\n",
            params.workerID, elapsed.hours, elapsed.minutes, elapsed.seconds);
        gaps_flush();
    }
    return result;
}

