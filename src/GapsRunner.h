#ifndef __COGAPS_GAPS_RUNNER_H__
#define __COGAPS_GAPS_RUNNER_H__

#include "GapsStatistics.h"
#include "GibbsSampler.h"

#include "data_structures/Matrix.h"

class GapsRunner
{
private:

    unsigned mNumRows;
    unsigned mNumCols;

    AmplitudeGibbsSampler mASampler;
    PatternGibbsSampler mPSampler;
    GapsStatistics mStatistics;

    bool mSamplePhase;

    unsigned mNumUpdatesA;
    unsigned mNumUpdatesP;

    void updateSampler(unsigned nA, unsigned nP, unsigned nCores);

public:

    GapsRunner(const RowMatrix &data, unsigned nPatterns, float alphaA,
        float alphaP, float maxGibbsMassA, float maxGibbsMassP,
        bool singleCell);

    GapsRunner(const std::string &pathToData, unsigned nPatterns, float alphaA,
        float alphaP, float maxGibbsMassA, float maxGibbsMassP,
        bool singleCell);
    
    void run(unsigned nIter, unsigned outputFreq, bool printMessages,
        unsigned nCores);

    unsigned nRow() const { return mNumRows; }
    unsigned nCol() const { return mNumCols; }

    void setUncertainty(const RowMatrix &S);
    void setUncertainty(const std::string &pathToMatrix);

    void startSampling() { mSamplePhase = true; }

    void displayStatus(unsigned outFreq, unsigned current, unsigned total);

    void setFixedMatrix(char which, const RowMatrix &mat);

    ColMatrix AMean() const { return mStatistics.AMean(); }
    RowMatrix PMean() const { return mStatistics.PMean(); }
    ColMatrix AStd() const { return mStatistics.AStd(); }
    RowMatrix PStd() const { return mStatistics.PStd(); }
};

#endif // __COGAPS_GAPS_RUNNER_H__