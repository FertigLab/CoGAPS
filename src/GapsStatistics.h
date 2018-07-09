#ifndef __COGAPS_GAPS_STATISTICS_H__
#define __COGAPS_GAPS_STATISTICS_H__

#include "GibbsSampler.h"
#include "data_structures/Matrix.h"

enum PumpThreshold
{
    PUMP_UNIQUE=1,
    PUMP_CUT=2
};

class GapsStatistics
{
private:

    ColMatrix mAMeanMatrix;
    ColMatrix mAStdMatrix;
    RowMatrix mPMeanMatrix;
    RowMatrix mPStdMatrix;
    
    unsigned mStatUpdates;
    unsigned mNumPatterns;

    ColMatrix mPumpMatrix;
    PumpThreshold mPumpThreshold;
    unsigned mPumpStatUpdates;

public:

    GapsStatistics(unsigned nRow, unsigned nCol, unsigned nPatterns,
        PumpThreshold t=PUMP_CUT);

    void update(const AmplitudeGibbsSampler &ASampler,
        const PatternGibbsSampler &PSampler);

    ColMatrix Amean() const;
    RowMatrix Pmean() const;
    ColMatrix Asd() const;
    RowMatrix Psd() const;

    float meanChiSq(const AmplitudeGibbsSampler &ASampler) const;

    // PUMP statistics
    void updatePump(const AmplitudeGibbsSampler &ASampler,
        const PatternGibbsSampler &PSampler);

    RowMatrix pumpMatrix() const;
    RowMatrix meanPattern();
    void patternMarkers(ColMatrix normedA, RowMatrix normedP, ColMatrix &statMatrix);

    // serialization
    friend Archive& operator<<(Archive &ar, GapsStatistics &stat);
    friend Archive& operator>>(Archive &ar, GapsStatistics &stat);
};

#endif