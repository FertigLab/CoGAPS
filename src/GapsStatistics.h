#ifndef __COGAPS_GAPS_STATISTICS_H__
#define __COGAPS_GAPS_STATISTICS_H__

#include "GibbsSampler.h"
#include "data_structures/Matrix.h"

class GapsStatistics
{
private:

    ColMatrix mAMeanMatrix;
    ColMatrix mAStdMatrix;
    ColMatrix mPMeanMatrix;
    ColMatrix mPStdMatrix;
    
    unsigned mStatUpdates;
    unsigned mNumPatterns;

public:

    GapsStatistics(unsigned nRow, unsigned nCol, unsigned nPatterns);

    void update(const GibbsSampler &ASampler, const GibbsSampler &PSampler);

    ColMatrix Amean() const;
    ColMatrix Pmean() const;
    ColMatrix Asd() const;
    ColMatrix Psd() const;

    float meanChiSq(const GibbsSampler &PSampler) const;

    // serialization
    friend Archive& operator<<(Archive &ar, GapsStatistics &stat);
    friend Archive& operator>>(Archive &ar, GapsStatistics &stat);
};

#endif