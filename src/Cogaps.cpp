#include "GibbsSampler.h"

#include <Rcpp.h>
#include <ctime>
#include <boost/date_time/posix_time/posix_time.hpp>

enum GapsPhase
{
    GAPS_CALIBRATION,
    GAPS_COOLING,
    GAPS_SAMPLING
};

static void runGibbsSampler(GibbsSampler &sampler, unsigned nIterTotal,
unsigned& nIterA, unsigned& nIterP, Vector& chi2Vec, Vector& aAtomVec,
Vector& pAtomVec, GapsPhase phase);

// [[Rcpp::export]]
Rcpp::List cogaps(Rcpp::NumericMatrix DMatrix, Rcpp::NumericMatrix SMatrix,
unsigned nFactor, double alphaA, double alphaP, unsigned nEquil,
unsigned nEquilCool, unsigned nSample, double maxGibbsMassA,
double maxGibbsMassP, int seed=-1, bool messages=false, bool singleCellRNASeq=false)
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

    // create the gibbs sampler
    GibbsSampler sampler(DMatrix, SMatrix, nFactor, alphaA, alphaP,
        maxGibbsMassA, maxGibbsMassP, singleCellRNASeq);

    // initial number of iterations for each matrix
    unsigned nIterA = 10;
    unsigned nIterP = 10;

    // calibration phase
    Vector chi2Vec(nEquil);
    Vector nAtomsAEquil(nEquil);
    Vector nAtomsPEquil(nEquil);
    runGibbsSampler(sampler, nEquil, nIterA, nIterP, chi2Vec,
        nAtomsAEquil, nAtomsPEquil, GAPS_CALIBRATION);

    // cooling phase
    Vector trash(nEquilCool);
    runGibbsSampler(sampler, nEquilCool, nIterA, nIterP, trash,
        trash, trash, GAPS_COOLING);

    // sampling phase
    Vector chi2VecSample(nSample);
    Vector nAtomsASample(nSample);
    Vector nAtomsPSample(nSample);
    runGibbsSampler(sampler, nSample, nIterA, nIterP, chi2VecSample,
        nAtomsASample, nAtomsPSample, GAPS_SAMPLING);

    // debug checks
#ifdef GAPS_DEBUG
    sampler.checkAtomicMatrixConsistency();
    sampler.checkAPMatrix();
#endif

    // combine chi2 vectors
    chi2Vec.concat(chi2VecSample);

    //Just leave the snapshots as empty lists
    return Rcpp::List::create(
        Rcpp::Named("Amean") = sampler.AMeanRMatrix(),
        Rcpp::Named("Asd") = sampler.AStdRMatrix(),
        Rcpp::Named("Pmean") = sampler.PMeanRMatrix(),
        Rcpp::Named("Psd") = sampler.PStdRMatrix(),
        Rcpp::Named("ASnapshots") = Rcpp::List::create(),
        Rcpp::Named("PSnapshots") = Rcpp::List::create(),
        Rcpp::Named("atomsAEquil") = nAtomsAEquil.rVec(),
        Rcpp::Named("atomsASamp") = nAtomsASample.rVec(),
        Rcpp::Named("atomsPEquil") = nAtomsPEquil.rVec(),
        Rcpp::Named("atomsPSamp") = nAtomsPSample.rVec(),
        Rcpp::Named("chiSqValues") = chi2Vec.rVec(),
        Rcpp::Named("randSeed") = seedUsed
#ifdef GAPS_DEBUG
        , Rcpp::Named("randTypes") = Rcpp::wrap(gaps::random::getTypes()),
        Rcpp::Named("randValues") = Rcpp::wrap(gaps::random::getValues()),
        Rcpp::Named("propHistoryA") = Rcpp::wrap(sampler.proposalHistory('A')),
        Rcpp::Named("propHistoryP") = Rcpp::wrap(sampler.proposalHistory('P')),
        Rcpp::Named("atomHistoryA") = Rcpp::wrap(sampler.atomHistory('A')),
        Rcpp::Named("atomHistoryP") = Rcpp::wrap(sampler.atomHistory('P'))
#endif
    );
}

static void runGibbsSampler(GibbsSampler &sampler, unsigned nIterTotal,
unsigned& nIterA, unsigned& nIterP, Vector& chi2Vec, Vector& aAtomVec,
Vector& pAtomVec, GapsPhase phase)
{
    double tempDenom = (double)nIterTotal / 2.0;

    for (unsigned i = 0; i < nIterTotal; ++i)
    {
        if (phase == GAPS_CALIBRATION)
        {
            sampler.setAnnealingTemp(std::min(((double)i + 1.0) / tempDenom, 1.0));
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
        }

        chi2Vec(i) = sampler.chi2();
        double numAtomsA = sampler.totalNumAtoms('A');
        double numAtomsP = sampler.totalNumAtoms('P');
        aAtomVec(i) = numAtomsA;
        pAtomVec(i) = numAtomsP;

        if (phase != GAPS_COOLING)
        {
            nIterA = gaps::random::poisson(std::max(numAtomsA, 10.0));
            nIterP = gaps::random::poisson(std::max(numAtomsP, 10.0));
        }
    }
}