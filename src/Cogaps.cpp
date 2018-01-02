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

static std::vector<Rcpp::NumericMatrix> snapshotsA;
static std::vector<Rcpp::NumericMatrix> snapshotsP;

static void runGibbsSampler(GibbsSampler &sampler, unsigned nIterTotal,
unsigned& nIterA, unsigned& nIterP, Vector& chi2Vec, Vector& aAtomVec,
Vector& pAtomVec, GapsPhase phase, unsigned numOutputs, bool messages,
unsigned numSnapshots);

// [[Rcpp::export]]
Rcpp::List cogaps(Rcpp::NumericMatrix DMatrix, Rcpp::NumericMatrix SMatrix,
unsigned nFactor, float alphaA, float alphaP, unsigned nEquil,
unsigned nEquilCool, unsigned nSample, float maxGibbsMassA,
float maxGibbsMassP, Rcpp::NumericMatrix fixedPatterns,
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

    // create the gibbs sampler
    GibbsSampler sampler(DMatrix, SMatrix, nFactor, alphaA, alphaP,
        maxGibbsMassA, maxGibbsMassP, singleCellRNASeq, fixedPatterns,
        whichMatrixFixed);

    // initial number of iterations for each matrix
    unsigned nIterA = 10;
    unsigned nIterP = 10;

    // calibration phase
    Vector chi2Vec(nEquil);
    Vector nAtomsAEquil(nEquil);
    Vector nAtomsPEquil(nEquil);
    runGibbsSampler(sampler, nEquil, nIterA, nIterP, chi2Vec,
        nAtomsAEquil, nAtomsPEquil, GAPS_CALIBRATION, numOutputs, messages,
        numSnapshots);

    // cooling phase
    Vector trash(nEquilCool);
    runGibbsSampler(sampler, nEquilCool, nIterA, nIterP, trash,
        trash, trash, GAPS_COOLING, numOutputs, messages, numSnapshots);

    // sampling phase
    Vector chi2VecSample(nSample);
    Vector nAtomsASample(nSample);
    Vector nAtomsPSample(nSample);
    runGibbsSampler(sampler, nSample, nIterA, nIterP, chi2VecSample,
        nAtomsASample, nAtomsPSample, GAPS_SAMPLING, numOutputs, messages,
        numSnapshots);

    // combine chi2 vectors
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
    float tempDenom = (float)nIterTotal / 2.f;

    for (unsigned i = 0; i < nIterTotal; ++i)
    {
        if (phase == GAPS_CALIBRATION)
        {
            sampler.setAnnealingTemp(std::min(((float)i + 2.f) / tempDenom, 1.f));
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

        float numAtomsA = sampler.totalNumAtoms('A');
        float numAtomsP = sampler.totalNumAtoms('P');
        aAtomVec[i] = numAtomsA;
        pAtomVec[i] = numAtomsP;
        chi2Vec[i] = sampler.chi2();

        if (phase != GAPS_COOLING)
        {
            if ((i + 1) % numOutputs == 0 && messages)
            {
                std::string temp = phase == GAPS_CALIBRATION ? "Equil: " : "Samp: ";
                std::cout << temp << i + 1 << " of " << nIterTotal << ", Atoms:"
                    << numAtomsA << "(" << numAtomsP << ") Chi2 = "
                    << sampler.chi2() << '\n';
            }

            nIterA = gaps::random::poisson(std::max(numAtomsA, 10.f));
            nIterP = gaps::random::poisson(std::max(numAtomsP, 10.f));
        }
    }
}