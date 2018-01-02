#ifndef __COGAPS_GIBBSSAMPLER_H__
#define __COGAPS_GIBBSSAMPLER_H__

#include "AtomicSupport.h"
#include "Matrix.h"

#include <vector>

class GibbsSampler
{
private:
#ifdef GAPS_INTERNAL_TESTS
public:
#endif

    AtomicSupport mADomain, mPDomain;

    TwoWayMatrix mDMatrix, mSMatrix, mAPMatrix;
    RowMatrix mPMatrix;
    ColMatrix mAMatrix;

    RowMatrix mPMeanMatrix, mPStdMatrix;
    ColMatrix mAMeanMatrix, mAStdMatrix;
    unsigned mStatUpdates;

    float mMaxGibbsMassA;
    float mMaxGibbsMassP;

    float mAnnealingTemp;
    float mChi2;

    bool mSingleCellRNASeq;

    unsigned mNumFixedPatterns;
    char mFixedMat;

    bool death(AtomicSupport &domain, AtomicProposal &proposal);
    bool birth(AtomicSupport &domain, AtomicProposal &proposal);
    bool move(AtomicSupport &domain, AtomicProposal &proposal);
    bool exchange(AtomicSupport &domain, AtomicProposal &proposal);

    bool evaluateChange(AtomicSupport &domain, const AtomicProposal &proposal,
        float threshold, bool accept=false);

    float computeDeltaLL(const MatrixChange &change);

    float getGibbsMass(const MatrixChange &change);

    void updateAPMatrix_A(unsigned row, unsigned col, float delta);
    void updateAPMatrix_P(unsigned row, unsigned col, float delta);
    void updateAPMatrix(const MatrixChange &change);

    bool canUseGibbs(const MatrixChange &ch);
    void setChi2(float chi2);

public:

    GibbsSampler(Rcpp::NumericMatrix D, Rcpp::NumericMatrix S, unsigned nFactor,
        float alphaA, float alphaP, float maxGibbsMassA, float maxGibbsMassP,
        bool singleCellRNASeq, Rcpp::NumericMatrix fixedPat, char whichMat);

    void update(char matrixLabel);

    uint64_t totalNumAtoms(char matrixLabel) const;
    void setAnnealingTemp(float temp);
    float chi2() const;

    Rcpp::NumericMatrix AMeanRMatrix() const;
    Rcpp::NumericMatrix AStdRMatrix() const;
    Rcpp::NumericMatrix PMeanRMatrix() const;
    Rcpp::NumericMatrix PStdRMatrix() const;

    void updateStatistics();

    Rcpp::NumericMatrix getNormedMatrix(char mat);
};

#endif
