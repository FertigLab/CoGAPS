#ifndef __COGAPS_DISPATCHER_H__
#define __COGAPS_DISPATCHER_H__

#include "GapsRunner.h"

#include <string>

struct GapsReturn
{
    ColMatrix Amean;
    ColMatrix Asd;
    RowMatrix Pmean;
    RowMatrix Psd;
    
    float meanChiSq;
    uint32_t seed;
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
    bool mSingleCellRnaSeq;

    std:vector<GapsRunner*> mRunners;

    void runOneCycle(unsigned k);

public:

    GapsDispatcher() : mNumPatterns(3), mMaxIterations(1000),
        mOutputFrequency(250), mSeed(0), mAlphaA(0.01), mAlphaP(0.01),
        mMaxGibbsMassA(100.f), mMaxGibbsMassP(100.f), mPrintMessages(true),
        mSingleCellRnaSeq(false)
    {}

    void setNumPatterns(unsigned n)     { mNumPatterns = n; }
    void setMaxIterations(unsigned n)   { mMaxIterations = n; }
    void setOutputFrequency(unsigned n) { mOutputFrequency = n; }
    void setSeed(unsigned seed)         { mSeed = seed; }

    void printMessages(bool print) { mPrintMessages = print; }
    void singleCellRNASeq(bool sc) { mSingleCellRnaSeq = sc; }

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

    void useDefaultUncertainty();
    void setUncertainty(const std::string &pathToMatrix);
    void setUncertainty(const RowMatrix &S);

    void loadData(const RowMatrix &D);
    void loadData(const std::string &pathToData);    
    
    GapsReturn run();
};

#endif