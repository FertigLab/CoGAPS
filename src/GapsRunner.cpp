#include "GapsAssert.h"
#include "GapsPrint.h"
#include "GapsRunner.h"
#include "math/Random.h"
#include "math/SIMD.h"

#ifdef __GAPS_R_BUILD__
#include <Rcpp.h>
#endif

void GapsRunner::printMessages(bool print)
{
    mPrintMessages = print;
}

void GapsRunner::setOutputFrequency(unsigned n)
{
    mOutputFrequency = n;
}

void GapsRunner::setSparsity(float alphaA, float alphaP, bool singleCell)
{
    mASampler.setSparsity(alphaA, singleCell);
    mPSampler.setSparsity(alphaP, singleCell);
}

void GapsRunner::setMaxGibbsMass(float maxA, float maxP)
{
    mASampler.setMaxGibbsMass(maxA);
    mPSampler.setMaxGibbsMass(maxP);
}

void GapsRunner::setFixedMatrix(char which, const RowMatrix &mat)
{
    mFixedMatrix = which;
    if (which == 'A')
    {
        mASampler.setMatrix(ColMatrix(mat));
        mASampler.recalculateAPMatrix();
        mPSampler.sync(mASampler);
    }
    else if (which == 'P')
    {
        mPSampler.setMatrix(mat);
        mPSampler.recalculateAPMatrix();
        mASampler.sync(mPSampler);
    }
}

void GapsRunner::startSampling()
{
    mSamplePhase = true;
}

void GapsRunner::run(unsigned nIter, unsigned nCores)
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

        if (mPrintMessages)
        {
            displayStatus(i, nIter);
        }

        if (mSamplePhase)
        {
            mStatistics.update(mASampler, mPSampler);
        }
    }
}

unsigned GapsRunner::nRow() const
{
    return mASampler.dataRows();
}

unsigned GapsRunner::nCol() const
{
    return mASampler.dataCols();
}

ColMatrix GapsRunner::Amean() const
{
    return mStatistics.Amean();
}

RowMatrix GapsRunner::Pmean() const
{
    return mStatistics.Pmean();
}

ColMatrix GapsRunner::Asd() const
{
    return mStatistics.Asd();
}

RowMatrix GapsRunner::Psd() const
{
    return mStatistics.Psd();
}

float GapsRunner::meanChiSq() const
{
    return mStatistics.meanChiSq(mASampler);
}

Archive& operator<<(Archive &ar, GapsRunner &runner)
{
    ar << runner.mASampler << runner.mPSampler << runner.mStatistics <<
        runner.mPrintMessages << runner.mOutputFrequency <<
        runner.mFixedMatrix << runner.mSamplePhase << runner.mNumUpdatesA <<
        runner.mNumUpdatesP;
    return ar;
}

Archive& operator>>(Archive &ar, GapsRunner &runner)
{
    ar >> runner.mASampler >> runner.mPSampler >> runner.mStatistics >>
        runner.mPrintMessages >> runner.mOutputFrequency >>
        runner.mFixedMatrix >> runner.mSamplePhase >> runner.mNumUpdatesA >>
        runner.mNumUpdatesP;
    return ar;
}

void GapsRunner::updateSampler(unsigned nA, unsigned nP, unsigned nCores)
{
    if (mFixedMatrix != 'A')
    {
        mNumUpdatesA += nA;
        mASampler.update(nA, nCores);
        if (mFixedMatrix != 'P')
        {
            mPSampler.sync(mASampler);
        }
    }

    if (mFixedMatrix != 'P')
    {
        mNumUpdatesP += nP;
        mPSampler.update(nP, nCores);
        if (mFixedMatrix != 'A')
        {
            mASampler.sync(mPSampler);
        }
    }
}

void GapsRunner::displayStatus(unsigned current, unsigned total)
{
    if (mOutputFrequency > 0 && ((current + 1) % mOutputFrequency) == 0)
    {
        gaps_printf("%d of %d, Atoms:%lu(%lu) Chi2 = %.2f\n", current + 1,
            total, mASampler.nAtoms(), mPSampler.nAtoms(), mASampler.chi2());
    }
}
