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

    GapsResult(unsigned nrow, unsigned ncol) : Amean(nrow, ncol),
        Asd(nrow, ncol), Pmean(nrow, ncol), Psd(nrow, ncol), meanChiSq(0.f),
        seed(0)
    {}

    void writeCsv(const std::string &path);
};

// should be agnostic to external caller (R/Python/CLI)
class GapsDispatcher
{
private:

    unsigned mNumPatterns;
    unsigned mMaxIterations;
    unsigned mOutputFrequency;
    uint32_t mSeed;

    float mAlphaA;
    float mAlphaP;
    float mMaxGibbsMassA;
    float mMaxGibbsMassP;

    bool mPrintMessages;
    bool mSingleCell;

    bool mDataIsLoaded;    

    unsigned mNumCoresPerSet;

    std::vector<GapsRunner*> mRunners;

    char mFixedMatrix;

    void runOneCycle(unsigned k);

    GapsDispatcher(const GapsDispatcher &p); // don't allow copies
    GapsDispatcher& operator=(const GapsDispatcher &p); // don't allow copies

public:

    explicit GapsDispatcher(uint32_t seed=0) : mNumPatterns(3),
        mMaxIterations(1000), mOutputFrequency(250), mSeed(seed), mAlphaA(0.01),
        mAlphaP(0.01), mMaxGibbsMassA(100.f), mMaxGibbsMassP(100.f),
        mPrintMessages(true), mSingleCell(false), mDataIsLoaded(false),
        mNumCoresPerSet(1), mFixedMatrix('N')
    {
        gaps::random::setSeed(mSeed);
    }

    ~GapsDispatcher()
    {
        for (unsigned i = 0; i < mRunners.size(); ++i)
        {
            delete mRunners[i];
        }
    }

    void setNumPatterns(unsigned n)     { mNumPatterns = n; }
    void setMaxIterations(unsigned n)   { mMaxIterations = n; }
    void setOutputFrequency(unsigned n) { mOutputFrequency = n; }
    void setNumCoresPerSet(unsigned n)  { mNumCoresPerSet = n; }

    void printMessages(bool print) { mPrintMessages = print; }
    void singleCell(bool sc) { mSingleCell = sc; }

    void setAlpha(float alphaA, float alphaP)
    {
        mAlphaA = alphaA;
        mAlphaP = alphaP;
    }

    void setMaxGibbsMass(float maxA, float maxP)
    {
        mMaxGibbsMassA = maxA;
        mMaxGibbsMassP = maxP;
    }

    void setFixedMatrix(char which, const RowMatrix &mat)
    {
        mFixedMatrix = which;
        mRunners[0]->setFixedMatrix(which, mat);
    }
    
    void loadCheckpointFile(const std::string &pathToCptFile);

    void setUncertainty(const RowMatrix &S);
    void setUncertainty(const std::string &pathToMatrix);

    void loadData(const RowMatrix &D);
    void loadData(const std::string &pathToData);
    
    GapsResult run();
};

#endif // __COGAPS_GAPS_DISPATCHER_H__