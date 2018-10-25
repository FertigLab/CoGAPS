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
    else
    {
        return runCoGAPSAlgorithm<DenseGibbsSampler>(data, params,
            uncertainty, randState);
    }
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
const Sampler &ASampler, const Sampler &PSampler, bpt::ptime startTime,
char phase, unsigned iter)
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
char phase, unsigned iter)
{
    if (params.printMessages && params.outputFrequency > 0
    && ((iter + 1) % params.outputFrequency) == 0)
    {
        bpt::time_duration diff = bpt_now() - startTime;
        double perComplete = estimatedPercentComplete(params, ASampler,
            PSampler, startTime, phase, iter);
        double nSecondsCurrent = diff.total_seconds();
        double nSecondsTotal = nSecondsCurrent / perComplete;

        unsigned elapsedSeconds = static_cast<unsigned>(nSecondsCurrent);
        unsigned totalSeconds = static_cast<unsigned>(nSecondsTotal);

        unsigned elapsedHours = elapsedSeconds / 3600;
        elapsedSeconds -= elapsedHours * 3600;
        unsigned elapsedMinutes = elapsedSeconds / 60;
        elapsedSeconds -= elapsedMinutes * 60;

        unsigned totalHours = totalSeconds / 3600;
        totalSeconds -= totalHours * 3600;
        unsigned totalMinutes = totalSeconds / 60;
        totalSeconds -= totalMinutes * 60;

        gaps_printf("%d of %d, Atoms: %lu(%lu), ChiSq: %.0f, Time: %02d:%02d:%02d / %02d:%02d:%02d\n",
            iter + 1, params.nIterations, ASampler.nAtoms(),
            PSampler.nAtoms(), PSampler.chiSq(), elapsedHours, elapsedMinutes,
            elapsedSeconds, totalHours, totalMinutes, totalSeconds);
        gaps_flush();
    }
}

template <class Sampler>
static void updateSampler(const GapsParameters &params, Sampler &ASampler,
Sampler &PSampler, unsigned nA, unsigned nP)
{
    if (!params.useFixedMatrix || params.whichFixedMatrix != 'A')
    {
        ASampler.update(nA, params.maxThreads);
        if (!params.useFixedMatrix || params.whichFixedMatrix != 'P')
        {
            PSampler.sync(ASampler, params.maxThreads);
        }
    }

    if (!params.useFixedMatrix || params.whichFixedMatrix != 'P')
    {
        PSampler.update(nP, params.maxThreads);
        if (!params.useFixedMatrix || params.whichFixedMatrix != 'A')
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
        }
        displayStatus(params, ASampler, PSampler, startTime, phase, currentIter);
    }
}

template <class Sampler>
static void processFixedMatrix(const GapsParameters &params, Sampler &ASampler,
Sampler &PSampler)
{
    // check if we're fixing a matrix
    if (params.useFixedMatrix)
    {
        switch (params.whichFixedMatrix)
        {
            GAPS_ASSERT(params.fixedMatrix.nCol() == params.nPatterns);
            case 'A' :
                GAPS_ASSERT(params.fixedMatrix.nRow() == params.nGenes);
                ASampler.setMatrix(params.fixedMatrix);
                break;
            case 'P' :
                GAPS_ASSERT(params.fixedMatrix.nRow() == params.nSamples);
                PSampler.setMatrix(params.fixedMatrix);
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
    return result;
}

