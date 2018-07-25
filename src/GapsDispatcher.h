#ifndef __COGAPS_GAPS_DISPATCHER_H__
#define __COGAPS_GAPS_DISPATCHER_H__

#include "GapsRunner.h"
#include "math/SIMD.h"

#include <string>

struct GapsResult
{
    ColMatrix Amean;
    ColMatrix Asd;
    RowMatrix Pmean;
    RowMatrix Psd;
    
    float meanChiSq;
    uint32_t seed;

    GapsResult(unsigned nrow, unsigned ncol, uint32_t rngSeed) :
        Amean(nrow, ncol), Asd(nrow, ncol), Pmean(nrow, ncol), Psd(nrow, ncol),
        meanChiSq(0.f), seed(rngSeed)
    {}

    void writeCsv(const std::string &path);
    void writeTsv(const std::string &path);
    void writeGct(const std::string &path);
};

// should be agnostic to external caller (R/Python/CLI)
class GapsDispatcher
{
private:

    uint32_t mSeed;
    unsigned mNumPatterns;

    unsigned mMaxIterations;
    unsigned mNumCoresPerSet;
    bool mPrintMessages;

    unsigned mCheckpointsCreated;
    char mPhase; // 'C' for calibration, 'S' for sample

    unsigned mCheckpointInterval;
    std::string mCheckpointOutFile;

    bool mInitialized;    
    std::vector<GapsRunner*> mRunners;

    void runOneCycle(unsigned k);
    void createCheckpoint() const;

    GapsDispatcher(const GapsDispatcher &p); // don't allow copies
    GapsDispatcher& operator=(const GapsDispatcher &p); // don't allow copies

    template <class DataType>
    void loadData(const DataType &data, bool transposeData);

    template <class DataType>
    void loadData(const DataType &data, bool transposeData, bool partitionRows,
        const std::vector<unsigned> &indices);

public:

    GapsDispatcher();
    ~GapsDispatcher();

    template <class DataType>
    void initialize(const DataType &data, bool transposeData,
        unsigned nPatterns, uint32_t seed);

    template <class DataType>
    void initialize(const DataType &data, bool transposeData,
        unsigned nPatterns, uint32_t seed, bool partitionRows,
        const std::vector<unsigned> &indices);

    template <class DataType>
    void initialize(const DataType &data, bool transposeData,
        const std::string &cptFile);

    template <class DataType>
    void setUncertainty(const DataType &unc, bool transposeData);

    template <class DataType>
    void setUncertainty(const DataType &unc, bool transposeData,
        bool partitionRows, const std::vector<unsigned> &indices);

    void setMaxIterations(unsigned n);
    void printMessages(bool print);
    void setOutputFrequency(unsigned n);
    void setSparsity(float alphaA, float alphaP, bool singleCell);
    void setMaxGibbsMass(float maxA, float maxP);

    void setAMatrix(const RowMatrix &mat);
    void setPMatrix(const RowMatrix &mat);
    void setFixedMatrix(char which);
    
    void setNumCoresPerSet(unsigned n);
    void setCheckpointInterval(unsigned n);
    void setCheckpointOutFile(const std::string &path);
    
    GapsResult run();
};

template <class DataType>
void GapsDispatcher::initialize(const DataType &data, bool transposeData,
unsigned nPatterns, uint32_t seed)
{
    mSeed = seed;
    mNumPatterns = nPatterns;
    gaps::random::setSeed(mSeed);
    
    loadData(data, transposeData);
    mInitialized = true;
}

template <class DataType>
void GapsDispatcher::initialize(const DataType &data, bool transposeData,
unsigned nPatterns, uint32_t seed, bool partitionRows,
const std::vector<unsigned> &indices)
{
    mSeed = seed;
    mNumPatterns = nPatterns;
    gaps::random::setSeed(mSeed);
    
    loadData(data, transposeData, partitionRows, indices);
    mInitialized = true;
}

template <class DataType>
void GapsDispatcher::initialize(const DataType &data, bool transposeData,
const std::string &cptFile)
{
    Archive ar(cptFile, ARCHIVE_READ);
    gaps::random::load(ar);

    ar >> mSeed >> mNumPatterns >> mMaxIterations >> mPrintMessages >>
        mCheckpointsCreated >> mPhase;

    loadData(data, transposeData);
    ar >> *mRunners[0];
    mInitialized = true;
}

template <class DataType>
void GapsDispatcher::setUncertainty(const DataType &unc, bool transposeData)
{
    GAPS_ASSERT(mInitialized);
    mRunners[0]->setUncertainty(unc, transposeData);
}

template <class DataType>
void GapsDispatcher::setUncertainty(const DataType &unc, bool transposeData,
bool partitionRows, const std::vector<unsigned> &indices)
{
    GAPS_ASSERT(mInitialized);
    mRunners[0]->setUncertainty(unc, transposeData, partitionRows, indices);
}

template <class DataType>
void GapsDispatcher::loadData(const DataType &data, bool transposeData)
{
    gaps_printf("Loading Data...");
    gaps_flush();
    mRunners.push_back(new GapsRunner(data, transposeData, mNumPatterns));
    gaps_printf("Done!\n");
}

template <class DataType>
void GapsDispatcher::loadData(const DataType &data, bool transposeData,
bool partitionRows, const std::vector<unsigned> &indices)
{
    gaps_printf("Loading Data...");
    gaps_flush();
    mRunners.push_back(new GapsRunner(data, transposeData, partitionRows,
        indices, mNumPatterns));
    gaps_printf("Done!\n");
}

#endif // __COGAPS_GAPS_DISPATCHER_H__