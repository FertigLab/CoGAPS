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

    double mMaxGibbsMassA;
    double mMaxGibbsMassP;

    double mAnnealingTemp;
    double mChi2;

    bool mSingleCellRNASeq;

    unsigned mNumFixedPatterns;
    char mFixedMat;

    bool death(AtomicSupport &domain, AtomicProposal &proposal);
    bool birth(AtomicSupport &domain, AtomicProposal &proposal);
    bool move(AtomicSupport &domain, AtomicProposal &proposal);
    bool exchange(AtomicSupport &domain, AtomicProposal &proposal);

    bool evaluateChange(AtomicSupport &domain, const AtomicProposal &proposal,
        double threshold, bool accept=false);

    double computeDeltaLL(const MatrixChange &change);

    double getGibbsMass(const MatrixChange &change);

    void updateAPMatrix_A(unsigned row, unsigned col, double delta);
    void updateAPMatrix_P(unsigned row, unsigned col, double delta);
    void updateAPMatrix(const MatrixChange &change);

    bool canUseGibbs(const MatrixChange &ch);
    void setChi2(double chi2);

public:

    GibbsSampler(Rcpp::NumericMatrix D, Rcpp::NumericMatrix S, unsigned nFactor,
        double alphaA, double alphaP, double maxGibbsMassA, double maxGibbsMassP,
        bool singleCellRNASeq, Rcpp::NumericMatrix fixedPat, char whichMat);

    void update(char matrixLabel);

    uint64_t totalNumAtoms(char matrixLabel) const;
    void setAnnealingTemp(double temp);
    double chi2() const;

    Rcpp::NumericMatrix AMeanRMatrix() const;
    Rcpp::NumericMatrix AStdRMatrix() const;
    Rcpp::NumericMatrix PMeanRMatrix() const;
    Rcpp::NumericMatrix PStdRMatrix() const;

    void updateStatistics();

    Rcpp::NumericMatrix getNormedMatrix(char mat);
};

#endif
