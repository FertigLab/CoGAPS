#ifndef __COGAPS_GAPS_RESULT__
#define __COGAPS_GAPS_RESULT__

#include "GapsStatistics.h"
#include "data_structures/Matrix.h"

#include <stdint.h>
#include <string>

struct GapsResult
{
    Matrix Amean;
    Matrix Asd;
    Matrix Pmean;
    Matrix Psd;
    Matrix pumpMatrix;
    Matrix meanPatternAssignment;

    std::vector<float> chisqHistory;
    std::vector<unsigned> atomHistoryA;
    std::vector<unsigned> atomHistoryP;

    uint64_t totalUpdates;
    uint32_t seed;
    unsigned totalRunningTime;

    float meanChiSq;
    float averageQueueLengthA;
    float averageQueueLengthP;

    explicit GapsResult(const GapsStatistics &stat);

    void writeToFile(const std::string &fullPath);
    void writeCsv(const std::string &path);
    void writeTsv(const std::string &path);
    void writeGct(const std::string &path);
};

#endif // __COGAPS_GAPS_RESULT__