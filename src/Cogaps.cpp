#include "GibbsSampler.h"

#include <Rcpp.h>
#include <ctime>
#include <boost/date_time/posix_time/posix_time.hpp>

typedef std::vector<Rcpp::NumericMatrix> SnapshotList;

enum GapsPhase
{
    GAPS_CALIBRATION,
    GAPS_COOLING,
    GAPS_SAMPLING
};

// declaration order of member variables important!
// initialization depends on it
struct GapsInternalState
{
    GibbsSampler sampler;
    
    Vector chi2VecEquil;
    Vector nAtomsAEquil;
    Vector nAtomsPEquil;

    Vector chi2VecSample;
    Vector nAtomsASample;
    Vector nAtomsPSample;

    unsigned nIterA;
    unsigned nIterP;
    
    unsigned nEquil;
    unsigned nEquilCool;
    unsigned nSample;

    unsigned iteration;
    GapsPhase phase;

    SnapshotList snapshotsA;
    SnapshotList snapshotsP;

    GapsInternalState(Rcpp::NumericMatrix DMatrix, Rcpp::NumericMatrix SMatrix,
        unsigned nFactor, double alphaA, double alphaP, unsigned nE,
        unsigned nEC, unsigned nS, double maxGibbsMassA,
        double maxGibbsMassP, Rcpp::NumericMatrix fixedPatterns,
        char whichMatrixFixed, bool messages, bool singleCellRNASeq,
        unsigned numOutputs, unsigned numSnapshots)
            :
        chi2VecEquil(nEquil), nAtomsAEquil(nEquil), nAtomsPEquil(nEquil),
        chi2VecSample(nSample), nAtomsASample(nSample), nAtomsPSample(nSample),
        nIterA(10), nIterP(10), nEquil(nE), nEquilCool(nEC), nSample(nS),
        iteration(0), phase(GAPS_CALIBRATION),
        sampler(DMatrix, SMatrix, nFactor, alphaA, alphaP,
            maxGibbsMassA, maxGibbsMassP, singleCellRNASeq, fixedPatterns,
            whichMatrixFixed)
    {}

    GibbsSampler ;

    GapsInternalState(std::ifstream &file)
        :
    {}


};

// [[Rcpp::export]]
Rcpp::List cogapsFromCheckpoint(const std::string &fileName)
{   
    // open file
    std::ifstream file(fileName);

    // seed random number generator

    // load state from file and run
    GapsInternalState state(file);
    return runCogaps(state);
}

// [[Rcpp::export]]
Rcpp::List cogaps(Rcpp::NumericMatrix DMatrix, Rcpp::NumericMatrix SMatrix,
unsigned nFactor, double alphaA, double alphaP, unsigned nEquil,
unsigned nEquilCool, unsigned nSample, double maxGibbsMassA,
double maxGibbsMassP, Rcpp::NumericMatrix fixedPatterns,
char whichMatrixFixed, int seed, bool messages, bool singleCellRNASeq,
unsigned numOutputs, unsigned numSnapshots)
{
    // set seed
    uint32_t seedUsed = static_cast<uint32_t>(seed);
    if (seed < 0)
    {
        boost::posix_time::ptime epoch(boost::gregorian::date(1970,1,1)); 
        boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
        boost::posix_time::time_duration diff = now - epoch;
        seedUsed = static_cast<uint32_t>(diff.total_milliseconds() % 1000);
    }
    gaps::random::setSeed(seedUsed);

    GapsInternalState state(DMatrix, SMatrix, nFactor, alphaA, alphaP,
        nEquil, nEquilCool, nSample, maxGibbsMassA, maxGibbsMassP,
        fixedPatterns, whichMatrixFixed, messages, singleCellRNASeq,
        numOutputs, numSnapshots);
    return runCogaps(state);
}

Rcpp::List runCogaps(GapsInternalState &state)
{
    // calibration phase
    runGibbsSampler(sampler, nEquil, nIterA, nIterP, chi2Vec,
        nAtomsAEquil, nAtomsPEquil, GAPS_CALIBRATION, numOutputs, messages,
        numSnapshots);

    // cooling phase
    Vector trash(nEquilCool);
    runGibbsSampler(sampler, nEquilCool, nIterA, nIterP, trash,
        trash, trash, GAPS_COOLING, numOutputs, messages, numSnapshots);

    // sampling phase
    runGibbsSampler(sampler, nSample, nIterA, nIterP, chi2VecSample,
        nAtomsASample, nAtomsPSample, GAPS_SAMPLING, numOutputs, messages,
        numSnapshots);

    // combine chi2 vectors
    Vector chi2Vec(chi2VecEquil);
    chi2Vec.concat(chi2VecSample);

    return Rcpp::List::create(
        Rcpp::Named("Amean") = sampler.AMeanRMatrix(),
        Rcpp::Named("Asd") = sampler.AStdRMatrix(),
        Rcpp::Named("Pmean") = sampler.PMeanRMatrix(),
        Rcpp::Named("Psd") = sampler.PStdRMatrix(),
        Rcpp::Named("ASnapshots") = Rcpp::wrap(snapshotsA),
        Rcpp::Named("PSnapshots") = Rcpp::wrap(snapshotsP),
        Rcpp::Named("atomsAEquil") = nAtomsAEquil.rVec(),
        Rcpp::Named("atomsASamp") = nAtomsASample.rVec(),
        Rcpp::Named("atomsPEquil") = nAtomsPEquil.rVec(),
        Rcpp::Named("atomsPSamp") = nAtomsPSample.rVec(),
        Rcpp::Named("chiSqValues") = chi2Vec.rVec(),
        Rcpp::Named("randSeed") = seedUsed
    );
}

static void runGibbsSampler(GibbsSampler &sampler, unsigned nIterTotal,
unsigned& nIterA, unsigned& nIterP, Vector& chi2Vec, Vector& aAtomVec,
Vector& pAtomVec, GapsPhase phase, unsigned numOutputs, bool messages,
unsigned numSnapshots)
{
    double tempDenom = (double)nIterTotal / 2.0;

    for (unsigned i = 0; i < nIterTotal; ++i)
    {
        // TODO check if checkpoint should be created

        if (phase == GAPS_CALIBRATION)
        {
            sampler.setAnnealingTemp(std::min(((double)i + 2.0) / tempDenom, 1.0));
        }

        for (unsigned j = 0; j < nIterA; ++j)
        {
            sampler.update('A');
        }

        for (unsigned j = 0; j < nIterP; ++j)
        {
            sampler.update('P');
        }

        if (phase == GAPS_SAMPLING)
        {
            sampler.updateStatistics();
            if (numSnapshots > 0 && (i + 1) % (nIterTotal / numSnapshots) == 0)
            {
                snapshotsA.push_back(sampler.getNormedMatrix('A'));
                snapshotsP.push_back(sampler.getNormedMatrix('P'));
            }
        }

        double numAtomsA = sampler.totalNumAtoms('A');
        double numAtomsP = sampler.totalNumAtoms('P');
        aAtomVec(i) = numAtomsA;
        pAtomVec(i) = numAtomsP;
        chi2Vec(i) = sampler.chi2();

        if (phase != GAPS_COOLING)
        {
            if ((i + 1) % numOutputs == 0 && messages)
            {
                std::string temp = phase == GAPS_CALIBRATION ? "Equil: " : "Samp: ";
                std::cout << temp << i + 1 << " of " << nIterTotal << ", Atoms:"
                    << numAtomsA << "(" << numAtomsP << ") Chi2 = "
                    << sampler.chi2() << '\n';
            }

            nIterA = gaps::random::poisson(std::max(numAtomsA, 10.0));
            nIterP = gaps::random::poisson(std::max(numAtomsP, 10.0));
        }
    }
}

