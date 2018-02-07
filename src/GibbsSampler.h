#ifndef __COGAPS_GIBBSSAMPLER_H__
#define __COGAPS_GIBBSSAMPLER_H__

#include "AtomicSupport.h"
#include "Matrix.h"

#include <vector>

enum PumpThreshold
{
    PUMP_UNIQUE=1,
    PUMP_CUT=2
};

class GibbsSampler
{
private:
#ifdef GAPS_INTERNAL_TESTS
public:
#endif

    TwoWayMatrix mDMatrix, mSMatrix, mAPMatrix;

    ColMatrix mAMatrix;
    RowMatrix mPMatrix;

    AtomicSupport mADomain, mPDomain;

    ColMatrix mAMeanMatrix, mAStdMatrix;
    RowMatrix mPMeanMatrix, mPStdMatrix;
    unsigned mStatUpdates;

    ColMatrix mPumpMatrix;
    PumpThreshold mPumpThreshold;
    unsigned mPumpStatUpdates;

    float mMaxGibbsMassA;
    float mMaxGibbsMassP;

    float mAnnealingTemp;

    bool mSingleCellRNASeq;

    unsigned mNumFixedPatterns;
    char mFixedMat;

    void death(AtomicSupport &domain, AtomicProposal &proposal);
    void birth(AtomicSupport &domain, AtomicProposal &proposal);
    void move(AtomicSupport &domain, AtomicProposal &prop);
    void exchange(AtomicSupport &domain, AtomicProposal &proposal);

    void evaluateChange(AtomicSupport &domain, const AtomicProposal &proposal,
        MatrixChange &change, float threshold, bool accept=false);

    float computeDeltaLL(const MatrixChange &change);

    float getGibbsMass(const MatrixChange &change);

    void updateAPMatrix_A(unsigned row, unsigned col, float delta);
    void updateAPMatrix_P(unsigned row, unsigned col, float delta);
    void updateAPMatrix(const MatrixChange &change);

    bool canUseGibbs(const MatrixChange &ch);

public:

    GibbsSampler(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
        unsigned nFactor);
    GibbsSampler(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
        unsigned nFactor, float alphaA, float alphaP, float maxGibbmassA,
        float maxGibbmassP, bool singleCellRNASeq, char whichMatrixFixed,
        const Rcpp::NumericMatrix &FP, PumpThreshold pumpThreshold);

    void update(char matrixLabel);

    uint64_t totalNumAtoms(char matrixLabel) const;
    void setAnnealingTemp(float temp);
    float chi2() const;

    ColMatrix normedAMatrix() const;
    RowMatrix normedPMatrix() const;

    unsigned nRow() const {return mDMatrix.nRow();}
    unsigned nCol() const {return mDMatrix.nCol();}
    unsigned nFactor() const {return mAMatrix.nCol();}

    // statistics
    void updateStatistics();
    void updatePumpStatistics();
    Rcpp::NumericMatrix AMeanRMatrix() const;
    Rcpp::NumericMatrix AStdRMatrix() const;
    Rcpp::NumericMatrix PMeanRMatrix() const;
    Rcpp::NumericMatrix PStdRMatrix() const;
    Rcpp::NumericMatrix pumpMatrix() const;
    Rcpp::NumericMatrix meanPattern();
    void patternMarkers(RowMatrix normedA, RowMatrix normedP, ColMatrix &statMatrix);
    float meanChiSq() const;

    // serialization
    friend Archive& operator<<(Archive &ar, GibbsSampler &sampler);
    friend Archive& operator>>(Archive &ar, GibbsSampler &sampler);
};

#endif
