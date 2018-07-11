#ifndef __COGAPS_GAPS_RUNNER_H__
#define __COGAPS_GAPS_RUNNER_H__

#include "GapsStatistics.h"
#include "GibbsSampler.h"

#include "data_structures/Matrix.h"

// boost time helpers
#include <boost/date_time/posix_time/posix_time.hpp>
namespace bpt = boost::posix_time;
#define bpt_now() bpt::microsec_clock::local_time()

class GapsRunner
{
private:
    
    AmplitudeGibbsSampler mASampler;
    PatternGibbsSampler mPSampler;
    GapsStatistics mStatistics;

    bool mPrintMessages;
    unsigned mOutputFrequency;
    char mFixedMatrix;
    bool mSamplePhase;

    unsigned mNumUpdatesA;
    unsigned mNumUpdatesP;
        
    bpt::ptime mStartTime;

    void updateSampler(unsigned nA, unsigned nP, unsigned nCores);
    void displayStatus(unsigned current, unsigned total);

    double estimatePercentComplete();

public:

    template <class DataType>
    GapsRunner(const DataType &data, unsigned nPatterns);

    template <class DataType>
    void setUncertainty(const DataType &unc);

    void printMessages(bool print);
    void setOutputFrequency(unsigned n);
    void setSparsity(float alphaA, float alphaP, bool singleCell);
    void setMaxGibbsMass(float maxA, float maxP);

    void setFixedMatrix(char which, const RowMatrix &mat);

    void startSampling();

    void run(unsigned nIter, unsigned nCores);

    unsigned nRow() const;
    unsigned nCol() const;

    ColMatrix Amean() const;
    RowMatrix Pmean() const;
    ColMatrix Asd() const;
    RowMatrix Psd() const;
    float meanChiSq() const;

    void startClock();

    // serialization
    friend Archive& operator<<(Archive &ar, GapsRunner &runner);
    friend Archive& operator>>(Archive &ar, GapsRunner &runner);
};

template <class DataType>
GapsRunner::GapsRunner(const DataType &data, unsigned nPatterns)
    :
mASampler(data, nPatterns), mPSampler(data, nPatterns),
mStatistics(mASampler.dataRows(), mPSampler.dataCols(), nPatterns),
mSamplePhase(false), mNumUpdatesA(0), mNumUpdatesP(0), mFixedMatrix('N')
{
    mASampler.sync(mPSampler);
    mPSampler.sync(mASampler);
}

template <class DataType>
void GapsRunner::setUncertainty(const DataType &unc)
{
    mASampler.setUncertainty(unc);
    mPSampler.setUncertainty(unc);
}

#endif // __COGAPS_GAPS_RUNNER_H__