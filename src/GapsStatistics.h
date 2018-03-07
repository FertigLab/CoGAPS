#ifndef __GAPS_GAPS_STATISTICS_H__
#define __GAPS_GAPS_STATISTICS_H__

#include "GibbsSampler.h"
#include "Matrix.h"

class GapsStatistics
{
private:

    ColMatrix mAMeanMatrix;
    ColMatrix mAStdMatrix;
    RowMatrix mPMeanMatrix;
    RowMatrix mPStdMatrix;
    
    unsigned mStatUpdates;    

public:

    GapsStatistics(unsigned nRow, unsigned nCol, unsigned nFactor);

    Rcpp::NumericMatrix AMean();
    Rcpp::NumericMatrix AStd();
    Rcpp::NumericMatrix PMean();
    Rcpp::NumericMatrix PStd();

    void update(const AmplitudeGibbsSampler &ASampler,
        const PatternGibbsSampler &PSampler);
};

#endif