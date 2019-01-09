#include "GapsRunner.h"

#include "utils/Archive.h"

#ifdef __GAPS_R_BUILD__
#include <Rcpp.h>
#endif

#ifdef __GAPS_OPENMP__
#include <omp.h>
#endif

// boost time helpers
#include <boost/date_time/posix_time/posix_time.hpp>
namespace bpt = boost::posix_time;
#define bpt_now() bpt::microsec_clock::local_time()

// for conditionally printing status messages
#define GAPS_MESSAGE(b, m) \
    do { \
        if (b) { \
            gaps_printf(m); \
        } \
    } while(0)

struct GapsTime
{
    unsigned hours;
    unsigned minutes;
    unsigned seconds;

    GapsTime(unsigned totalSeconds)
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

// helper function, this dispatches the correct run function depending
// on the type of GibbsSampler needed for the given parameters
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

    if (params.useSparseOptimization)
    {
        return runCoGAPSAlgorithm<SparseGibbsSampler>(data, params,
            uncertainty, randState);
    }
    return runCoGAPSAlgorithm<DenseGibbsSampler>(data, params,
        uncertainty, randState);
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
const Sampler &ASampler, const Sampler &PSampler, char phase, unsigned iter)
{
    double nIter = static_cast<double>(iter);
    double nAtomsA = static_cast<double>(ASampler.nAtoms());
    double nAtomsP = static_cast<double>(PSampler.nAtoms());
    
    if (phase == 'S')
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
char phase, unsigned iter, GapsStatistics &stats)
{
    if (params.printMessages && params.outputFrequency > 0
    && ((iter + 1) % params.outputFrequency) == 0)
    {
        bpt::time_duration diff = bpt_now() - startTime;
        double perComplete = estimatedPercentComplete(params, ASampler,
            PSampler, phase, iter);

        GapsTime elapsedTime(static_cast<unsigned>(diff.total_seconds()));
        GapsTime totalTime(static_cast<unsigned>(diff.total_seconds() / perComplete));

        float cs = PSampler.chiSq();
        unsigned nA = ASampler.nAtoms();
        unsigned nP = PSampler.nAtoms();

        stats.addChiSq(cs);
        stats.addAtomCount(nA, nP);

        gaps_printf("%d of %d, Atoms: %d(%d), ChiSq: %.0f, Time: %02d:%02d:%02d / %02d:%02d:%02d\n",
            iter + 1, params.nIterations, nA, nP, cs, elapsedTime.hours,
            elapsedTime.minutes, elapsedTime.seconds, totalTime.hours,
            totalTime.minutes, totalTime.seconds);
        gaps_flush();
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
const GapsStatistics &stats, const GapsRng &rng, char phase, unsigned iter)
{
    if (params.checkpointInterval > 0
    && ((iter + 1) % params.checkpointInterval) == 0
    && !params.subsetData)
    {
        // create backup file
        std::rename(params.checkpointOutFile.c_str(),
            (params.checkpointOutFile + ".backup").c_str());
    
        // create checkpoint file
        Archive ar(params.checkpointOutFile, ARCHIVE_WRITE);
        ar << params;
        ar << *randState;
        ar << ASampler << PSampler << stats << phase << iter << rng;
        
        // delete backup file
        std::remove((params.checkpointOutFile + ".backup").c_str());

        ASampler.extraInitialization();
        PSampler.extraInitialization();
    }
}

template <class Sampler>
static void processCheckpoint(GapsParameters &params, Sampler &ASampler,
Sampler &PSampler, GapsRandomState *randState, GapsStatistics &stats,
GapsRng &rng, char &phase, unsigned &currentIter)
{    
    // check if running from checkpoint, get all saved data
    if (params.useCheckPoint)
    {
        Archive ar(params.checkpointFile, ARCHIVE_READ);
        ar >> params;
        ar >> *randState;
        ar >> ASampler >> PSampler >> stats >> phase >> currentIter >> rng;
    }
}

template <class Sampler>
static void runOnePhase(const GapsParameters &params, Sampler &ASampler,
Sampler &PSampler, GapsStatistics &stats, const GapsRandomState *randState,
GapsRng &rng, bpt::ptime startTime, char phase, unsigned &currentIter)
{
    for (; currentIter < params.nIterations; ++currentIter)
    {
        #ifdef __GAPS_R_BUILD__
        Rcpp::checkUserInterrupt();
        #endif

        createCheckpoint(params, ASampler, PSampler, randState, stats,
            rng, phase, currentIter);
        
        // set annealing temperature in calibration phase
        if (phase == 'C')
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

        if (phase == 'S')
        {
            stats.update(ASampler, PSampler);
            if (params.takePumpSamples)
            {
                stats.updatePump(ASampler);
            }
        }
        displayStatus(params, ASampler, PSampler, startTime, phase,
            currentIter, stats);
    }
}

template <class Sampler>
static void processFixedMatrix(const GapsParameters &params, Sampler &ASampler,
Sampler &PSampler)
{
    // check if we're fixing a matrix
    if (params.useFixedPatterns)
    {
        switch (params.whichMatrixFixed)
        {
            GAPS_ASSERT(params.fixedPatterns.nCol() == params.nPatterns);
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
    GAPS_MESSAGE(params.printMessages, "Running in debug mode\n");
    #endif

    // load data into gibbs samplers
    // we transpose the data in the A sampler so that the update step
    // is symmetrical for each sampler, this simplifies the code 
    // within the sampler, note the subsetting genes/samples flag must be
    // flipped if we are flipping the transpose flag
    GAPS_MESSAGE(params.printMessages, "Loading Data...");
    Sampler ASampler(data, !params.transposeData, !params.subsetGenes,
        params.alphaA, params.maxGibbsMassA, params, randState);
    Sampler PSampler(data, params.transposeData, params.subsetGenes,
        params.alphaP, params.maxGibbsMassP, params, randState);
    processUncertainty(params, ASampler, PSampler, uncertainty);
    processFixedMatrix(params, ASampler, PSampler);
    GAPS_MESSAGE(params.printMessages, "Done!\n");

    // these variables will get overwritten by checkpoint if provided
    GapsStatistics stats(params.nGenes, params.nSamples, params.nPatterns);
    GapsRng rng(randState);
    char phase = 'C';
    unsigned currentIter = 0;
    processCheckpoint(params, ASampler, PSampler, randState, stats, rng,
        phase, currentIter);
    calculateNumberOfThreads(params);

    // sync samplers and run any additional initialization needed
    ASampler.sync(PSampler);
    PSampler.sync(ASampler);
    ASampler.extraInitialization();
    PSampler.extraInitialization();

    // record start time
    bpt::ptime startTime = bpt_now();

    // fallthrough through phases, allows algorithm to be resumed in any phase
    GAPS_ASSERT(phase == 'C' || phase == 'S');
    switch (phase)
    {
        case 'C':
            GAPS_MESSAGE(params.printMessages, "-- Calibration Phase --\n");
            runOnePhase(params, ASampler, PSampler, stats, randState, rng,
                startTime, phase, currentIter);
            phase = 'S';
            currentIter = 0;

        case 'S':
            GAPS_MESSAGE(params.printMessages, "-- Sampling Phase --\n");
            runOnePhase(params, ASampler, PSampler, stats, randState, rng,
                startTime, phase, currentIter);
    }
    
    // get result
    GapsResult result(stats);
    result.meanChiSq = stats.meanChiSq(PSampler);
    result.averageQueueLengthA = ASampler.getAverageQueueLength();
    result.averageQueueLengthP = PSampler.getAverageQueueLength();

    if (params.takePumpSamples)
    {
        result.pumpMatrix = stats.pumpMatrix();
        result.meanPatternAssignment = stats.meanPattern();
    }

    // if we are running distributed, each worker needs to print when it's done
    if (params.runningDistributed)
    {
        bpt::time_duration diff = bpt_now() - startTime;
        GapsTime elapsed(static_cast<unsigned>(diff.total_seconds()));
        gaps_printf("    worker %d is finished! Time: %02d:%02d:%02d\n",
            params.workerID, elapsed.hours, elapsed.minutes, elapsed.seconds);
        gaps_flush();
    }

    return result;
}

