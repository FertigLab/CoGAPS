#ifndef __COGAPS_GAPS_RUNNER_H__
#define __COGAPS_GAPS_RUNNER_H__

#include "GapsStatistics.h"
#include "GibbsSampler.h"

#include "data_structures/Matrix.h"

// boost time helpers
#include <boost/date_time/posix_time/posix_time.hpp>
namespace bpt = boost::posix_time;
#define bpt_now() bpt::microsec_clock::local_time()

struct GapsResult
{
    ColMatrix Amean;
    ColMatrix Asd;
    RowMatrix Pmean;
    RowMatrix Psd;
    
    float meanChiSq;
    uint32_t seed;

    GapsResult(const GapsStatistics &stat) :
        Amean(stat.Amean()), Asd(stat.Asd()), Pmean(stat.Pmean()),
        Psd(stat.Psd()), meanChiSq(0.f), seed(0)
    {}

    void writeCsv(const std::string &path);
    void writeTsv(const std::string &path);
    void writeGct(const std::string &path);
};

class GapsRunner
{
private:
    
    AmplitudeGibbsSampler mASampler;
    PatternGibbsSampler mPSampler;
    GapsStatistics mStatistics;

    char mFixedMatrix;
    unsigned mMaxIterations;
    
    unsigned mMaxThreads;
    bool mPrintMessages;
    unsigned mOutputFrequency;
    std::string mCheckpointOutFile;
    unsigned mCheckpointInterval;

    bpt::ptime mStartTime;
    char mPhase;
    unsigned mCurrentIteration;

    // only kept since they need to be written to the start of every checkpoint
    unsigned mNumPatterns;
    uint32_t mSeed;

    unsigned mNumUpdatesA;
    unsigned mNumUpdatesP;
        
    void runOnePhase();
    void updateSampler(unsigned nA, unsigned nP);
    double estimatedPercentComplete() const;
    void displayStatus();
    void createCheckpoint();

public:

    template <class DataType>
    GapsRunner(const DataType &data, bool transposeData, unsigned nPatterns,
        bool partitionRows, const std::vector<unsigned> &indices);

    template <class DataType>
    void setUncertainty(const DataType &unc, bool transposeData,
        bool partitionRows, const std::vector<unsigned> &indices);

    void setFixedMatrix(char which, const Matrix &mat);

    void recordSeed(uint32_t seed);
    uint32_t getSeed() const;

    void setMaxIterations(unsigned nIterations);
    void setSparsity(float alphaA, float alphaP, bool singleCell);
    void setMaxGibbsMass(float maxA, float maxP);
    
    void setMaxThreads(unsigned nThreads);
    void setPrintMessages(bool print);
    void setOutputFrequency(unsigned n);
    void setCheckpointOutFile(const std::string &outFile);
    void setCheckpointInterval(unsigned interval);

    GapsResult run(bool printThreads=true);

    // serialization
    friend Archive& operator>>(Archive &ar, GapsRunner &runner);
};

// problem with passing file parser - need to read it twice
template <class DataType>
GapsRunner::GapsRunner(const DataType &data, bool transposeData,
unsigned nPatterns, bool partitionRows, const std::vector<unsigned> &indices)
    :
mASampler(data, transposeData, nPatterns, partitionRows, indices),
mPSampler(data, transposeData, nPatterns, partitionRows, indices),
mStatistics(mASampler.dataRows(), mASampler.dataCols(), nPatterns),
mFixedMatrix('N'), mMaxIterations(1000), mMaxThreads(1), mPrintMessages(true),
mOutputFrequency(500), mCheckpointOutFile("gaps_checkpoint.out"),
mCheckpointInterval(0), mPhase('C'), mCurrentIteration(0),
mNumPatterns(nPatterns), mSeed(0), mNumUpdatesA(0), mNumUpdatesP(0)
{
    mASampler.sync(mPSampler);
    mPSampler.sync(mASampler);
}

template <class DataType>
void GapsRunner::setUncertainty(const DataType &unc, bool transposeData,
bool partitionRows, const std::vector<unsigned> &indices)
{
    mASampler.setUncertainty(unc, transposeData, partitionRows, indices);
    mPSampler.setUncertainty(unc, transposeData, partitionRows, indices);
}

#endif // __COGAPS_GAPS_RUNNER_H__