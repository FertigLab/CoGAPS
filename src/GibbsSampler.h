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

    TwoWayMatrix mDMatrix, mSMatrix, mAPMatrix;

    ColMatrix mAMatrix;
    RowMatrix mPMatrix;

    AtomicSupport mADomain, mPDomain;

    ColMatrix mAMeanMatrix, mAStdMatrix;
    RowMatrix mPMeanMatrix, mPStdMatrix;

    ColMatrix mPumpMatrix;
    Vector mLP;

    unsigned mStatUpdates;

    float mMaxGibbsMassA;
    float mMaxGibbsMassP;

    float mAnnealingTemp;

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

public:

    GibbsSampler(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
        unsigned nFactor);
    GibbsSampler(const Rcpp::NumericMatrix &D,
        const Rcpp::NumericMatrix &S, unsigned nFactor, float alphaA,
        float alphaP, float maxGibbmassA, float maxGibbmassP,
        bool singleCellRNASeq, char whichMatrixFixed,
        const Rcpp::NumericMatrix &FP);

    void update(char matrixLabel);

    uint64_t totalNumAtoms(char matrixLabel) const;
    void setAnnealingTemp(float temp);
    float chi2() const;

    Rcpp::NumericMatrix AMeanRMatrix() const;
    Rcpp::NumericMatrix AStdRMatrix() const;
    Rcpp::NumericMatrix PMeanRMatrix() const;
    Rcpp::NumericMatrix PStdRMatrix() const;
    float meanChiSq() const;

    void updateStatistics();

    Rcpp::NumericMatrix getNormedMatrix(char mat);

    unsigned nRow() const {return mDMatrix.nRow();}
    unsigned nCol() const {return mDMatrix.nCol();}
    unsigned nFactor() const {return mAMatrix.nCol();}

    friend Archive& operator<<(Archive &ar, GibbsSampler &sampler);
    friend Archive& operator>>(Archive &ar, GibbsSampler &sampler);
};

#endif
