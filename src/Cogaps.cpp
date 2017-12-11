#include "GibbsSampler.h"

#include <Rcpp.h>
#include <ctime>

void runGibbsSampler(GibbsSampler &sampler, unsigned iterations,
    bool updateTemp=false);

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
    runGibbsSampler(sampler, nEquil, true);
    runGibbsSampler(sampler, nEquilCool);
    runGibbsSampler(sampler, nSample);

    //Just leave the snapshots as empty lists
    return Rcpp::List::create(Rcpp::Named("randSeed") = seedUsed,
                              Rcpp::Named("numAtomA") = sampler.totalNumAtoms('A'));
}

void runGibbsSampler(GibbsSampler &sampler, unsigned iterations, bool updateTemp)
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

        double tempChiSq = sampler.chi2();
        uint64_t numAtomsA = sampler.totalNumAtoms('A');
        uint64_t numAtomsP = sampler.totalNumAtoms('P');

        nIterA = gaps::random::poisson(std::max(numAtomsA, (uint64_t)10));
        nIterP = gaps::random::poisson(std::max(numAtomsP, (uint64_t)10));
    }

}