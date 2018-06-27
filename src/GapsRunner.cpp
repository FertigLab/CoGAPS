#include "GapsAssert.h"
#include "GapsPrint.h"
#include "GapsRunner.h"
#include "math/Random.h"
#include "math/SIMD.h"

#ifdef __GAPS_R_BUILD__
#include <Rcpp.h>
#endif

#include <algorithm>

GapsRunner::GapsRunner(const RowMatrix &data, unsigned nPatterns, float alphaA,
float alphaP, float maxGibbsMassA, float maxGibbsMassP, bool singleCell)
    :
mNumRows(data.nRow()), mNumCols(data.nCol()),
mASampler(data, nPatterns, alphaA, maxGibbsMassA, singleCell),
mPSampler(data, nPatterns, alphaP, maxGibbsMassP, singleCell),
mStatistics(data.nRow(), data.nCol(), nPatterns),
mSamplePhase(false), mNumUpdatesA(0), mNumUpdatesP(0)
{
    mASampler.sync(mPSampler);
    mPSampler.sync(mASampler);
}

GapsRunner::GapsRunner(const std::string &pathToData, unsigned nPatterns,
float alphaA, float alphaP, float maxGibbsMassA, float maxGibbsMassP,
bool singleCell)
    :
mNumRows(FileParser(pathToData).nRow()), mNumCols(FileParser(pathToData).nCol()),
mASampler(pathToData, nPatterns, alphaA, maxGibbsMassA, singleCell),
mPSampler(pathToData, nPatterns, alphaP, maxGibbsMassP, singleCell),
mStatistics(mNumRows, mNumCols, nPatterns),
mSamplePhase(false), mNumUpdatesA(0), mNumUpdatesP(0)
{
    mASampler.sync(mPSampler);
    mPSampler.sync(mASampler);
}

void GapsRunner::setFixedMatrix(char which, const RowMatrix &mat)
{
    if (which == 'A')
    {
        mASampler.setMatrix(ColMatrix(mat));
    }
    else if (which == 'P')
    {
        mPSampler.setMatrix(mat);
    }
}

void GapsRunner::setUncertainty(const std::string &pathToMatrix)
{
    mASampler.setUncertainty(pathToMatrix);
    mPSampler.setUncertainty(pathToMatrix);
}

void GapsRunner::displayStatus(unsigned outFreq, unsigned current, unsigned total)
{
    if (outFreq > 0 && ((current + 1) % outFreq) == 0)
    {
        gaps_printf("%d of %d, Atoms:%lu(%lu) Chi2 = %.2f\n", current + 1,
            total, mASampler.nAtoms(), mPSampler.nAtoms(), mASampler.chi2());
    }
}

void GapsRunner::run(unsigned nIter, unsigned outputFreq, bool printMessages,
unsigned nCores)
{
    for (unsigned i = 0; i < nIter; ++i)
    {
        #ifdef __GAPS_R_BUILD__
        Rcpp::checkUserInterrupt();
        #endif
        
        if (!mSamplePhase)
        {        
            float temp = static_cast<float>(2 * i) / static_cast<float>(nIter);
            mASampler.setAnnealingTemp(std::min(1.f, temp));
            mPSampler.setAnnealingTemp(std::min(1.f, temp));
        }
    
        // number of updates per iteration is poisson 
        unsigned nA = gaps::random::poisson(std::max(mASampler.nAtoms(), 10ul));
        unsigned nP = gaps::random::poisson(std::max(mPSampler.nAtoms(), 10ul));
        updateSampler(nA, nP, nCores);

        if (printMessages)
        {
            displayStatus(outputFreq, i, nIter);
        }

        if (mSamplePhase)
        {
            mStatistics.update(mASampler, mPSampler);
        }
    }
}

void GapsRunner::updateSampler(unsigned nA, unsigned nP, unsigned nCores)
{
    mNumUpdatesA += nA;
    mASampler.update(nA, nCores);
    mPSampler.sync(mASampler);

    mNumUpdatesP += nP;
    mPSampler.update(nP, nCores);
    mASampler.sync(mPSampler);
}

void GapsRunner::setUncertainty(const RowMatrix &S)
{
    mASampler.setUncertainty(S);
    mPSampler.setUncertainty(S);
}
