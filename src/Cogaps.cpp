#include "GibbsSampler.h"

#include <Rcpp.h>
#include <ctime>

static void runGibbsSampler(GibbsSampler &sampler, unsigned iterations,
Vector &chi2Vec, Vector &aAtomVec, Vector &pAtomVec, bool updateTemp=false,
bool updateStats=false);

// [[Rcpp::export]]
Rcpp::List cogaps(Rcpp::NumericMatrix DMatrix, Rcpp::NumericMatrix SMatrix,
unsigned nFactor, double alphaA, double alphaP, unsigned nEquil,
unsigned nEquilCool, unsigned nSample, double maxGibbsMassA,
double maxGibbsMassP, int seed=-1, bool messages=false, bool singleCellRNASeq=false)
{
    // set seed
    uint32_t seedUsed = seed >= 0 ? static_cast<uint32_t>(seed)
        : static_cast<uint32_t>(std::time(0));
    gaps::random::setSeed(seedUsed);

    // create the gibbs sampler
    GibbsSampler sampler(DMatrix, SMatrix, nFactor, alphaA, alphaP,
        maxGibbsMassA, maxGibbsMassP, singleCellRNASeq);

    // run the sampler for each phase of the algorithm

    Vector chi2Vec(nEquil);
    Vector nAtomsAEquil(nEquil);
    Vector nAtomsPEquil(nEquil);
    runGibbsSampler(sampler, nEquil, chi2Vec, nAtomsAEquil, nAtomsPEquil, true);

    Vector trash(nEquilCool);
    runGibbsSampler(sampler, nEquilCool, trash, trash, trash);

    Vector chi2VecSample(nSample);
    Vector nAtomsASample(nSample);
    Vector nAtomsPSample(nSample);
    runGibbsSampler(sampler, nSample, chi2VecSample, nAtomsASample, nAtomsPSample,
        false, true);

    chi2Vec.concat(chi2VecSample);

    // compute statistics

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
        Rcpp::Named("randSeed") = seedUsed);
}

static void runGibbsSampler(GibbsSampler &sampler, unsigned iterations,
Vector& chi2Vec, Vector& aAtomVec, Vector& pAtomVec, bool updateTemp,
bool updateStats)
{
    unsigned nIterA = 10;
    unsigned nIterP = 10;
        
    double tempDenom = (double)iterations / 2.0;

    for (unsigned i = 0; i < iterations; ++i)
    {
        if (updateTemp)
        {
            sampler.setAnnealingTemp(std::min((i + 1) / tempDenom, 1.0));
        }

        for (unsigned j = 0; j < nIterA; ++j)
        {
            sampler.update('A');
        }

        for (unsigned j = 0; j < nIterP; ++j)
        {
            sampler.update('P');
        }

        if (updateStats)
        {
            sampler.updateStatistics();
        }

        chi2Vec(i) = sampler.chi2();
        uint64_t numAtomsA = sampler.totalNumAtoms('A');
        uint64_t numAtomsP = sampler.totalNumAtoms('P');
        aAtomVec(i) = (double) numAtomsA;
        pAtomVec(i) = (double) numAtomsP;

        nIterA = gaps::random::poisson(std::max(numAtomsA, (uint64_t)10));
        nIterP = gaps::random::poisson(std::max(numAtomsP, (uint64_t)10));
    }
}