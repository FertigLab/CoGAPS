#ifndef __COGAPS_GAPS_RESULT__
#define __COGAPS_GAPS_RESULT__

#include "data_structures/Matrix.h"

#include <stdint.h>
#include <string>
#include <vector>

class GapsStatistics;

struct GapsResult
{
    explicit GapsResult(const GapsStatistics &stat);
    explicit GapsResult();
    void writeToFile(const std::string &path);

    Matrix Amean;
    Matrix Asd;
    Matrix Pmean;
    Matrix Psd;
    Matrix pumpMatrix;
    Matrix meanPatternAssignment;
    std::vector<Matrix> equilibrationSnapshotsA;
    std::vector<Matrix> equilibrationSnapshotsP;
    std::vector<Matrix> samplingSnapshotsA;
    std::vector<Matrix> samplingSnapshotsP;
    std::vector<float> chisqHistory;
    std::vector<unsigned> atomHistoryA;
    std::vector<unsigned> atomHistoryP;
    uint64_t totalUpdates;
    uint32_t seed;
    unsigned totalRunningTime;
    float meanChiSq;
    float averageQueueLengthA;
    float averageQueueLengthP;
};

#endif // __COGAPS_GAPS_RESULT__