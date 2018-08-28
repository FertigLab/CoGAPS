#ifndef __COGAPS_GAPS_RESULT__
#define __COGAPS_GAPS_RESULT__

#include "GapsStatistics.h"
#include "data_structures/Matrix.h"

#include <stdint.h>
#include <string>

struct GapsResult
{
    ColMatrix Amean;
    ColMatrix Asd;
    ColMatrix Pmean;
    ColMatrix Psd;
    
    float meanChiSq;
    uint32_t seed;

    GapsResult(const GapsStatistics &stat);

    void writeToFile(const std::string &fullPath);
    void writeCsv(const std::string &path);
    void writeTsv(const std::string &path);
    void writeGct(const std::string &path);
};

#endif // __COGAPS_GAPS_RESULT__