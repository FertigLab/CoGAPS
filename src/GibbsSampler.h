#ifndef __COGAPS_GIBBSSAMPLER_H__
#define __COGAPS_GIBBSSAMPLER_H__

#include "AtomicSupport.h"
#include "Matrix.h"

#include <vector>

class GibbsSampler
{
private:

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

#ifdef GAPS_INTERNAL_TESTS
public:
#endif

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

#ifndef GAPS_INTERNAL_TESTS
public:
#endif

    GibbsSampler(Rcpp::NumericMatrix D, Rcpp::NumericMatrix S, unsigned nFactor,
        double alphaA, double alphaP, double maxGibbsMassA, double maxGibbsMassP,
        bool singleCellRNASeq);

    void update(char matrixLabel);

    uint64_t totalNumAtoms(char matrixLabel) const;
    void setAnnealingTemp(double temp);
    double chi2() const;

    Rcpp::NumericMatrix AMeanRMatrix() const;
    Rcpp::NumericMatrix AStdRMatrix() const;
    Rcpp::NumericMatrix PMeanRMatrix() const;
    Rcpp::NumericMatrix PStdRMatrix() const;

    void updateStatistics();

#ifdef GAPS_DEBUG
    void checkAtomicMatrixConsistency() const;
    void checkAPMatrix() const;
    std::vector<char> proposalTypeHistory(char label)
    {
        return label == 'A' ? mADomain.proposalTypeHistory()
            : mPDomain.proposalTypeHistory();
    }

    std::vector<double> proposalDelta1History(char label)
    {
        return label == 'A' ? mADomain.proposalDelta1History()
            : mPDomain.proposalDelta1History();
    }
    std::vector<double> proposalDelta2History(char label)
    {
        return label == 'A' ? mADomain.proposalDelta2History()
            : mPDomain.proposalDelta2History();
    }
    std::vector<char> acceptTypeHistory(char label)
    {
        return label == 'A' ? mADomain.acceptTypeHistory()
            : mPDomain.acceptTypeHistory();
    }

    std::vector<double> acceptDelta1History(char label)
    {
        return label == 'A' ? mADomain.acceptDelta1History()
            : mPDomain.acceptDelta1History();
    }
    std::vector<double> acceptDelta2History(char label)
    {
        return label == 'A' ? mADomain.acceptDelta2History()
            : mPDomain.acceptDelta2History();
    }
    std::vector<unsigned> atomHistory(char label)
    {
        return label == 'A' ? mADomain.atomHistory()
            : mPDomain.atomHistory();
    }
#endif
};

#endif
